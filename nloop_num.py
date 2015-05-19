from snappy import *
from random import randrange
from collections import OrderedDict
import multiprocessing

class NeumannZagierDatum():
    """
    description of class, comment entire class
    """
    def __init__(self, manifold, engine=None, verbose=False, file_name=None):
        self.manifold = manifold
        self.engine = engine
        self.verbose = verbose
        self.file_name = file_name

        self._raw_gluing_equations = manifold.gluing_equations()

        self.num_shapes = self._raw_gluing_equations.ncols() // 3
        self.num_eqns = self._raw_gluing_equations.nrows()
        self._eliminated_shapes = self.num_shapes * [1, ]

        self.computed_ptolemy = False

        self.nz = None

        pari.set_real_precision(100)


    def every_third_elem(L, shift):
        return [L[i] for i in range(shift, len(L), 3)]


    def epsilon(CC, d):
        return 2**(-CC.precision()//d)


    def all_shape_parameters(self, z):
        return [z, 1 / (1 - z), 1 - 1 / z]


    def in_threes(self, L):
        return [L[3 * i : 3 * (i + 1)] for i in range(len(L) // 3)]


    def shift_in_threes(self, L, shifts):
        return sum([X[s:] + X[:s] for X, s in zip(self.in_threes(L),\
            shifts)], [])


    def is_geom(self, M, c):
        vol = M.volume()
        for v in c.volume_numerical():
            if abs(v - vol) < 1e-10:
                return True
        return False


    def gluing_equations(self):
        eqns = self._raw_gluing_equations
        new_cols = self.shift_in_threes(eqns.columns(), [(i - 1) % 3 for\
            i in self._eliminated_shapes])
        return matrix(new_cols).transpose()


    def ABCbar(self):
        eqns = self.gluing_equations()
        n = self.num_shapes
        return [eqns.matrix_from_columns(range(i, 3 * n, 3))\
            for i in range(3)]


    def target_vector(self):
        """
        Answer times pi*i is right-hand side of
        gluing equations.
        """
        m = self.num_eqns
        c = self.manifold.num_cusps()
        return vector(ZZ, [2 for i in range(m - 2 * c)] + (2 * c) * [0])


    def ABv(self):
        A, B, C = self.ABCbar()
        one = vector(B.base_ring(), B.ncols() * [1])
        return A - B, C - B, self.target_vector() - B * one


    def shapes(self, precision=None, shapes=None):
        if shapes:
            base_shapes = shapes
        else:
            if precision == None:
                base_shapes = [CC(z) for z in\
                    self.manifold.tetrahedra_shapes(part='rect')]
            else:
                base_shapes = hypertorsion.snap.\
                    polished_tetrahedra_shapes(self.manifold,\
                    precision)

        return [self.all_shape_parameters(z)[(i - 1) % 3]
                 for z, i in zip(base_shapes, self._eliminated_shapes)]


    def ABv_square(self, rows_to_eliminate=None):
        A, B, v = self.ABv()
        M = block_matrix([[B, A, v.column()]])
        c = self.manifold.num_cusps()

        rows = range(M.nrows())[:-2 * c]
        rows += [2 * i + rows[-1] + 1 for i in range(c)]

        M = M.matrix_from_rows(rows)
        M = M.hermite_form(include_zero_rows=False)
        n = A.ncols()
        return M.matrix_from_columns(range(n, 2 * n)),\
            M.matrix_from_columns(range(n)), M.columns()[-1]


    def f_and_fddot(self):
        A, B, v = self.ABv_square()
        n = A.ncols()
        M = block_matrix([[A, B]])
        S, U, V = M.smith_form()
        d = S.diagonal()
        f = V * vector(ZZ, [x / y for x, y in zip(U * v, d)] + n * [0])
        assert M * f == v
        return vector(f[:n]), vector(f[n:])


    def make_B_nondegenerate(self):
        while det(self.ABv_square()[1]) == 0:
            self._eliminated_shapes = [randrange(3) for\
                i in range(self.num_shapes)]


    def compute_ptolemy_field_and_embedding(self):
        vol = self.manifold.volume()
        p = self.manifold.ptolemy_variety(2, 'all')
        if self.computed_ptolemy == False:
            if self.engine == "retrieve":
                try:
                    s = p.retrieve_solutions(verbose=self.verbose)\
                        .flatten(depth=2)
                except:
                    s = p.compute_solutions(engine=self.engine,\
                        verbose=self.verbose).flatten(depth=2)
            else:
                s = p.compute_solutions(engine=self.engine,\
                    verbose=self.verbose).flatten(depth=2)
            self.computed_ptolemy = s
        else:
            s = self.computed_ptolemy
        for sol in s:
            rsol = zip(pari('polroots(%s)' %sol.number_field()),\
                sol.numerical())
        for root, numerical_sol in rsol:
            if abs(vol - numerical_sol.volume_numerical()) < 1e-10:
                return sol.number_field(), root


    def check(self):
        shapes = vector(self.all_log_shapes())
        CC = shapes[0].parent()
        eqns = self.gluing_equations()
        pi_I = CC.pi() * CC.gen()
        error = eqns * shapes - pi_I * self.target_vector()
        assert error.norm(Infinity) < epsilon(CC, 2)

        A, B, v = self.ABv()
        z = vector(every_third_elem(shapes, 0))
        z_ddot = vector(every_third_elem(shapes, 2))
        error = A * z + B * z_ddot - pi_I * v
        assert error.norm(Infinity) < epsilon(CC, 2)


    def generate_nz_data(self):
        self.make_B_nondegenerate()
        temp_ABv = self.ABv_square()
        A = temp_ABv[0]
        B = temp_ABv[1]
        nu = temp_ABv[2]
        temp_fs = self.f_and_fddot()
        f = temp_fs[0]
        f_ddot = temp_fs[1]
        new_shapes = self.shapes()
        pol, embedding = self.compute_ptolemy_field_and_embedding()
        new_nz = (A, B, nu, f, f_ddot, pol, new_shapes, embedding)
        self.nz = new_nz
        if self.file_name != None:
            save(new_nz, self.file_name)


class nloop():
    """
    Compute the n-loop invariant S_n.

    Reference: ``The Quantum Content of the Gluing Equations'' by
    Dimofte and Garoufalidis.
    """
    def __init__(self, nzdata, n, diagrams):
        """Initializes class variables."""
        (A, B, nu, f, f_ddot, _, zees, _) = nzdata
        self.A = A
        self.B = B
        self.nu = nu
        self.f = f
        self.f_ddot = f_ddot
        self.zees = zees
        self.CC = self.zees[0].parent()
        self.prec = self.zees[0].prec()
        self.n = n
        self.diagrams = [g for g in diagrams if \
            self.feynman_loop_number(g) <= self.n]
        self.ver_factor = None

        self.prev = OrderedDict()


    def exponentiate_list(self, L, E):
        return prod([l ** e for l, e in zip(L, E)])


    def one_loop(self, precision=None, shapes=None):
        CC = self.zees[0].parent()
        shapes_dd = [1 - 1 / z for z in self.zees]
        D1 = diagonal_matrix(shapes_dd)
        D2 = diagonal_matrix([1 / z for z in self.zees])
        return (1 / CC(2)) * det(self.A * D1 + self.B * D2) *\
            self.exponentiate_list(self.zees, self.f_ddot) *\
            self.exponentiate_list(shapes_dd, -self.f)


    def pre_comp_polylog(self, index, z):
        """
        This function contains the polylogs commonly used
        in the calculation of the n-loop invariant.
        These polylogs are used a large number of times
        in a computation so they are saved to minimize waste.
        """
        if index == 1:
            return -ln(1 - z)
        if index == 0:
            return z / (1 - z)
        if index == -1:
            return z / (z ** 2 - 2 * z + 1)
        if index == -2:
            return (-z ** 2 - z) / (z ** 3 - 3 * z ** 2 + 3 * z - 1)
        if index == -3:
            return (z ** 3 + 4 * z ** 2 + z) / (z ** 4 - 4 * z ** 3 +\
                6 * z ** 2 - 4 * z + 1)
        if index == -4:
            return (-z ** 4 - 11 * z ** 3 - 11 * z ** 2 - z) / (z **\
                5 - 5 * z ** 4 + 10 * z ** 3 - 10 * z ** 2 + 5 * z - 1)
        if index == -5:
            return (z ** 5 + 26 * z ** 4 + 66 * z ** 3 + 26 * z ** 2 +\
                z) / (z ** 6 - 6 * z ** 5 + 15 * z ** 4 - 20 * z ** 3 +\
                15 * z ** 2 - 6 * z + 1)
        if index == -6:
            return (-z ** 6 - 57 * z ** 5 - 302 * z ** 4 - 302 * z **\
                3 - 57 * z ** 2 - z) / (z ** 7 - 7 * z ** 6 + 21 *\
                z ** 5 - 35 * z ** 4 + 35 * z ** 3 -\
                21 * z ** 2 + 7 * z - 1)
        if index == -7:
            return (z ** 7 + 120 * z ** 6 + 1191 * z ** 5 + 2416 *\
                z ** 4 + 1191 * z ** 3 + 120 * z ** 2 + z) / (z **\
                8 - 8 * z ** 7 + 28 * z ** 6 - 56 * z ** 5 + 70 *\
                z ** 4 - 56 * z ** 3 + 28 * z ** 2 - 8 * z + 1)
        if index == -8:
            return (-z ** 8 - 247 * z ** 7 - 4293 * z ** 6 -\
                15619 * z ** 5 - 15619 * z ** 4 - 4293 * z **\
                3 - 247 * z ** 2 - z) / (z ** 9 - 9 * z ** 8 +\
                36 * z ** 7 - 84 * z ** 6 + 126 * z ** 5 - 126 *\
                z ** 4 + 84 * z ** 3 - 36 * z ** 2 + 9 * z - 1)
        if index == -9:
            return (z ** 9 + 502 * z ** 8 + 14608 * z ** 7 +\
                88234 * z ** 6 + 156190 * z ** 5 + 88234 * z ** 4 +\
                14608 * z ** 3 + 502 * z ** 2 + z) / (z ** 10 -\
                10 * z ** 9 + 45 * z ** 8 - 120 * z ** 7 + 210 * z **\
                6 - 252 * z ** 5 + 210 * z ** 4 - 120 * z ** 3 + 45 *\
                z ** 2 - 10 * z + 1)
        if index == -10:
            return (-z ** 10 - 1013 * z ** 9 - 47840 * z ** 8 -\
                455192 * z ** 7 - 1310354 * z ** 6 - 1310354 * z ** 5 -\
                455192 * z ** 4 - 47840 * z ** 3 - 1013 * z ** 2 - z) /\
                (z ** 11 - 11 * z ** 10 + 55 * z ** 9 - 165 * z ** 8 +\
                330 * z ** 7 - 462 * z ** 6 + 462 * z ** 5 - 330 * z **\
                4 + 165 * z ** 3 - 55 * z ** 2 + 11 * z - 1)
        if index == -11:
            return (z ** 11 + 2036 * z ** 10 + 152637 * z ** 9 +\
                2203488 * z ** 8 + 9738114 * z**7 + 15724248 *\
                z ** 6 + 9738114 * z ** 5 + 2203488 * z ** 4 +\
                152637 * z ** 3 + 2036 * z ** 2 + z) / (z ** 12 - 12 *\
                z ** 11 + 66 * z ** 10 - 220 * z ** 9 + 495 * z ** 8 -\
                792 * z ** 7 + 924 * z ** 6 - 792 * z ** 5 + 495 * z **\
                4 - 220 * z ** 3 + 66 * z ** 2 - 12 * z + 1)
        if index == -12:
            return (-z ** 12 - 4083 * z ** 11 - 478271 * z ** 10 -\
                10187685 * z ** 9 - 66318474 * z ** 8 - 162512286 *\
                z ** 7 - 162512286 * z ** 6 - 66318474 * z ** 5 -\
                10187685 * z ** 4 - 478271 * z ** 3 - 4083 * z **\
                2 - z) / (z ** 13 - 13 * z ** 12 + 78 * z ** 11 -\
                286 * z ** 10 + 715 * z ** 9 - 1287 * z ** 8 + 1716 *\
                z ** 7 - 1716 * z ** 6 + 1287 * z ** 5 - 715 * z **\
                4 + 286 * z ** 3 - 78 * z ** 2 + 13 * z - 1)
        if index == -13:
            return (z ** 13 + 8178 * z ** 12 + 1479726 * z ** 11 +\
                45533450 * z ** 10 + 423281535 * z ** 9 + 1505621508 *\
                z ** 8 + 2275172004 * z ** 7 + 1505621508 * z ** 6 +\
                423281535 * z ** 5 + 45533450 * z ** 4 + 1479726 *\
                z ** 3 + 8178 * z ** 2 + z) / (z ** 14 - 14 * z ** 13 +\
                91 * z ** 12 - 364 * z ** 11 + 1001 * z ** 10 - 2002 *\
                z ** 9 + 3003 * z ** 8 - 3432 * z ** 7 + 3003 * z **\
                6 - 2002 * z ** 5 + 1001 * z ** 4 - 364 * z ** 3 +\
                91 * z ** 2 - 14 * z + 1)
        if index == -14:
            return (-z ** 14 - 16369 * z ** 13 - 4537314 * z ** 12 -\
                198410786 * z ** 11 - 2571742175 * z ** 10 -\
                12843262863 * z ** 9 - 27971176092 * z ** 8 -\
                27971176092 * z ** 7 - 12843262863 * z ** 6 -\
                2571742175 * z ** 5 - 198410786 * z ** 4 - 4537314 *\
                z ** 3 - 16369 * z ** 2 - z) / (z ** 15 - 15 * z **\
                14 + 105 * z ** 13 - 455 * z ** 12 + 1365 * z ** 11 -\
                3003 * z ** 10 + 5005 * z ** 9 - 6435 * z ** 8 +\
                6435 * z ** 7 - 5005 * z ** 6 + 3003 * z ** 5 -\
                1365 * z ** 4 + 455 * z ** 3 - 105 * z ** 2 + 15 *\
                z - 1)
        if index == -15:
            return (z ** 15 + 32752 * z ** 14 + 13824739 * z ** 13 +\
                848090912 * z ** 12 + 15041229521 * z ** 11 +\
                102776998928 * z ** 10 + 311387598411 * z ** 9 +\
                447538817472 * z ** 8 + 311387598411 * z ** 7 +\
                102776998928 * z ** 6 + 15041229521 * z ** 5 +\
                848090912 * z ** 4 + 13824739 * z ** 3 + 32752 * z **\
                2 + z) / (z ** 16 - 16 * z ** 15 + 120 * z ** 14 -\
                560 * z ** 13 + 1820 * z ** 12 - 4368 * z ** 11 +\
                8008 * z ** 10 - 11440 * z ** 9 + 12870 * z ** 8 -\
                11440 * z ** 7 + 8008 * z ** 6 - 4368 * z ** 5 +\
                1820 * z ** 4 - 560 * z ** 3 + 120 * z ** 2 -\
                16 * z + 1)
        if index == -16:
            return (-z ** 16 - 65519 * z ** 15 - 41932745 * z ** 14 -\
                3572085255 * z ** 13 - 85383238549 * z ** 12 -\
                782115518299 * z ** 11 - 3207483178157 * z ** 10 -\
                6382798925475 * z ** 9 - 6382798925475 * z ** 8 -\
                3207483178157 * z ** 7 - 782115518299 * z ** 6 -\
                85383238549 * z ** 5 - 3572085255 * z ** 4 -\
                41932745 * z ** 3 - 65519 * z ** 2 - z) / (z **\
                17 - 17 * z ** 16 + 136 * z ** 15 - 680 * z **\
                14 + 2380 * z ** 13 - 6188 * z ** 12 + 12376 *\
                z ** 11 - 19448 * z ** 10 + 24310 * z ** 9 - 24310 *\
                z ** 8 + 19448 * z ** 7 - 12376 * z ** 6 + 6188 * z **\
                5 - 2380 * z ** 4 + 680 * z ** 3 - 136 * z ** 2 +\
                17 * z - 1)
        if index == -17:
            return (z ** 17 + 131054 * z ** 16 + 126781020 * z ** 15 +\
                14875399450 * z ** 14 + 473353301060 * z ** 13 +\
                5717291972382 * z ** 12 + 31055652948388 * z ** 11 +\
                83137223185370 * z ** 10 + 114890380658550 * z ** 9 +\
                83137223185370 * z ** 8 + 31055652948388 * z ** 7 +\
                5717291972382 * z ** 6 + 473353301060 * z ** 5 +\
                14875399450 * z ** 4 + 126781020 * z ** 3 + 131054 *\
                z ** 2 + z) / (z ** 18 - 18 * z ** 17 + 153 * z **\
                16 - 816 * z ** 15 + 3060 * z ** 14 - 8568 * z **\
                13 + 18564 * z ** 12 - 31824 * z ** 11 + 43758 *\
                z ** 10 - 48620 * z ** 9 + 43758 * z ** 8 - 31824 *\
                z ** 7 + 18564 * z ** 6 - 8568 * z ** 5 + 3060 *\
                z ** 4 - 816 * z ** 3 + 153 * z ** 2 - 18 * z + 1)
        if index == -18:
            return (-z ** 18 - 262125 * z ** 17 - 382439924 * z ** 16 -\
                61403313100 * z ** 15 - 2575022097600 * z ** 14 -\
                40457344748072 * z ** 13 - 285997074307300 * z ** 12 -\
                1006709967915228 * z ** 11 - 1865385657780650 * z **\
                10 - 1865385657780650 * z ** 9 - 1006709967915228 *\
                z ** 8 - 285997074307300 * z ** 7 - 40457344748072 *\
                z ** 6 - 2575022097600 * z ** 5 - 61403313100 * z **\
                4 - 382439924 * z ** 3 - 262125 * z ** 2 - z) / (z **\
                19 - 19 * z ** 18 + 171 * z ** 17 - 969 * z ** 16 +\
                3876 * z ** 15 - 11628 * z ** 14 + 27132 * z ** 13 -\
                50388 * z ** 12 + 75582 * z ** 11 - 92378 * z ** 10 +\
                92378 * z ** 9 - 75582 * z ** 8 + 50388 * z ** 7 -\
                27132 * z ** 6 + 11628 * z ** 5 - 3876 * z ** 4 +\
                969 * z ** 3 - 171 * z ** 2 + 19 * z - 1)
        if index == -19:
            return (z ** 19 + 524268 * z ** 18 + 1151775897 * z ** 17 +\
                251732291184 * z ** 16 + 13796160184500 * z ** 15 +\
                278794377854832 * z ** 14 + 2527925001876036 * z **\
                13 + 11485644635009424 * z ** 12 + 27862280567093358 *\
                z ** 11 + 37307713155613000 * z ** 10 +\
                27862280567093358 * z ** 9 + 11485644635009424 * z **\
                8 + 2527925001876036 * z ** 7 + 278794377854832 * z **\
                6 + 13796160184500 * z ** 5 + 251732291184 * z ** 4 +\
                1151775897 * z ** 3 + 524268 * z ** 2 + z) / (z **\
                20 - 20 * z ** 19 + 190 * z ** 18 - 1140 * z ** 17 +\
                4845 * z ** 16 - 15504 * z ** 15 + 38760 * z ** 14 -\
                77520 * z ** 13 + 125970 * z ** 12 - 167960 * z **\
                11 + 184756 * z ** 10 - 167960 * z ** 9 + 125970 *\
                z ** 8 - 77520 * z ** 7 + 38760 * z ** 6 - 15504 *\
                z ** 5 + 4845 * z ** 4 - 1140 * z ** 3 + 190 * z **\
                2 - 20 * z + 1)


    def feynman_loop_number(self, diagram):
        """
        Calculate the Feynman Loop Number of a Diagram.

        The Feynman Loop Number of a connected looped multigraph
        is the number of 1-vertices+2-vertices + the number of loops
        """
        if diagram.num_edges() == 0:
            return 0
        return diagram.degree().count(1) + diagram.degree().count(2) +\
            diagram.num_edges() - diagram.num_verts() + 1


    def symmetry_factor(self, diag):
        """
        Calculate the symmetry factor of a diagram.

        This is equal to the order of the group of vertex
        permutations preserving edges times k! for each
        k-multiedge times 2^number of loops
        """
        symfactor = diag.automorphism_group().cardinality()
        for foo in diag.vertices():
            for bar in range(foo, diag.num_verts()):
                conecs = diag.adjacency_matrix()[foo][bar]
                symfactor = kronecker_delta(foo, bar) * symfactor * 2 **\
                    conecs * factorial(conecs) + (1 - kronecker_delta(\
                    foo, bar)) * symfactor * factorial(conecs)
        return QQ(1) / symfactor


    def bernoulli_plus_half(self, m):
        """Return bernoulli number with convention B1=+1/2."""
        return bernoulli(m) * (-1) ** m


    def polylogarithm(self, index, z):
        """Give the nth polylogarithm evaluated at z."""
        _ = gp.set_precision(self.prec)
        if (index, z) in self.prev:
            return self.prev[(index, z)]
        if index <= 1 and index >= -19:
            tmp = self.CC(self.pre_comp_polylog(index, z))
            self.prev[(index, z)] = tmp
            if len(self.prev) > 1000:
                self.prev.popitem(last=False)
            return tmp
        return self.CC(gp.subst(gp(polylog(index, x)), x, z))


    def gamma(self, eye, kay, ell):
        """Return the gamma equation for vertex_factor_tensor."""
        if kay == 0:
            return sum([self.polylogarithm(2 - self.n, 1 / z) for z \
                in self.zees]) * self.bernoulli_plus_half(self.n) /\
                factorial(self.n) + kronecker_delta(self.n, 2) * \
                (self.f * self.B.inverse() * self.A *\
                self.f / 8)[0][0]
        return (-1) ** kay * sum([h ** (bar - 1) / factorial(bar) *
            self.bernoulli_plus_half(bar) * self.polylogarithm(\
            2 - bar - kay, 1 / self.zees[eye]) for bar in \
            range((kronecker_delta(kay, 1) + kronecker_delta(kay, 2)),\
            1 + (kronecker_delta(kay, 1) + kronecker_delta(kay, 2)) + \
            self.n - ell)]) - kronecker_delta(kay, 1) * self.CC(0.5) * \
            (self.B.inverse() * self.nu)[eye]


    def vertex_factor_tensor(self):
        """
        Generate vertex gamma as a tensor access values.
        Output is in the form
        vertexgamma[feynman_loop_number][vertex_degree][ith_shape_parameter]
        """
        var('h')
        return [[[self.gamma(eye, kay, ell) for eye in
            range(len(self.zees))] for kay in range(2 * self.n + 1)]\
            for ell in range(self.n + 1)]


    def diagram_contribution_to_nloop(self, diagram):
        """The diagram contribution to the n-loop invariant."""
        N = len(self.zees)
        hamil = -self.B.inverse() * self.A + diagonal_matrix([1 / (1 - z)\
            for z in self.zees])
        prop = h * hamil.inverse()
        temp_sum = self.CC(0)
        for foo in range(N ** diagram.num_verts()):
            indices = [floor(foo / (N ** bar)) % N for bar in
                range(diagram.num_verts())]
            temp_sum += prod([prop[indices[eee[0]]]
                [indices[eee[1]]] for eee in diagram.edges
                (labels=False)] +
                [self.ver_factor[self.feynman_loop_number(diagram)]
                [diagram.degree()[vee]][indices[vee]] for vee in
                diagram.vertices()]).expand()
        return self.symmetry_factor(diagram) *\
            temp_sum.coeff(h, self.n - 1)


    def nloop_invariant(self):
        """The Dimofte-Garoufalidis n-loop invariant."""
        self.ver_factor = self.vertex_factor_tensor()

        PROCESSES = multiprocessing.cpu_count()
        #print 'cpu_count() = %d\n' % multiprocessing.cpu_count()
        pool = multiprocessing.Pool(PROCESSES)
        collect_results = pool.map(self.diagram_contribution_to_nloop,\
            self.diagrams)
        loop_invar = sum(collect_results) + self.ver_factor[self.n][0][0]
        return loop_invar


def nloop_from_manifold(manifold, n, diagrams, engine=None, verbose=False,\
    file_name=None):
    D = NeumannZagierDatum(manifold, engine, verbose, file_name)
    D.generate_nz_data()

    E = nloop(D.nz, n, all_diagrams)

    if n == 1:
        return [E.one_loop(), D.nz]
    return [E.nloop_invariant(), D.nz]


def nloop_from_nzdatum(nz, n, diagrams):
    E = nloop(nz, n, all_diagrams)

    if n == 1:
        return [E.one_loop(), nz]
    return [E.nloop_invariant(), nz]
