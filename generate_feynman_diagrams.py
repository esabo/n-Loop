"""
Generate all Feynman diagrams contributing to nloop invariant

Reference: ``The Quantum Content of the Gluing Equations'' by
Dimofte and Garoufalidis.
"""
from sets import Set
import copy

def holes(diagram):
    """
    Calculates the rank of the free group of a connected graph
    Equal to #edges - #vertices + 1
    """
    return diagram.num_edges()-diagram.num_verts()+1

def num_small_degree(diagram):
    """
    Calculates the number of vertices with degree 1 or 2.
    
    This is meant to work around a bug in python where the
    degree of a looped vertex in an immutable graph is 
    incorrectly computed.
    """
    diagram=diagram.copy(immutable=False)
    degree_list=diagram.degree()
    return degree_list.count(1)+degree_list.count(2)

def feynman_loop_number(diagram):
    """
    Calculates the Feynman Loop Number of diagram.
    The Feynman Loop Number of a connected looped multigraph
    is the number 
    # 1-vertices + # 2-vertices + # holes
    """
    if diagram.num_edges()==0:
        return 0
    else:
        return num_small_degree(diagram)+holes(diagram)
    
def loop_children(diagram, n):
    '''
    Returns a set of children diagrams obtained from diagram by adding
    an edge. Excludes some diagrams with all descendents having
    loop number greater than n.
    '''
    diagram=diagram.copy(immutable=False)
    if feynman_loop_number(diagram) > n:
        verts=[dex for dex,foo in enumerate(diagram.degree()) if foo==2]
    elif feynman_loop_number(diagram) == n:
        verts=[dex for dex,foo in enumerate(diagram.degree()) if foo<=2]
    else:
        verts=diagram.vertices()
    add_edges=Set([(foo,bar) for foo in verts for bar in diagram.vertices()])
    children=Set([])
    for edge in add_edges:
        child=diagram.copy()
        child.add_edge(edge)
        child=child.canonical_label()
        child=child.copy(immutable=True)
        children.add(child)
    return children

def nloop_diagrams(n):
    '''
    Generate the all connected looped multigraphs
    with Feynmann Loop Number at most n.
    
    The algorithm calls the Sage tree generator to generate trees.
    Loops and multiedges are added iteratively.           
    '''
    will_loop=Set([])
    done_loop=Set([])
    contributes=Set([])
    start_trees=[bar for foo in range(1,2*n-1) for bar in graphs.trees(foo)]
    for diagram in start_trees:
        if num_small_degree(diagram)<=2*n:
            diagram.allow_loops(True)
            diagram.allow_multiple_edges(True)
            diagram=diagram.canonical_label()
            diagram_copy=diagram.copy(immutable=True)
            will_loop.add(diagram_copy)
            if feynman_loop_number(diagram)<= n and feynman_loop_number(diagram)>0:
                contributes.add(diagram_copy)
    while len(will_loop)!=0:
#        print len(will_loop)
        diagram=will_loop.pop()
        children=loop_children(diagram,n)
        for child in children:
            if holes(child)+.5*num_small_degree(child)>n:
                done_loop.add(child)                
            if feynman_loop_number(child)<=n:
                contributes.add(child)
        done_loop.add(diagram)    
        children.difference_update(done_loop)
        will_loop.union_update(children)
    contributing_list=[diagram.copy(immutable=False) for diagram in contributes]
    contributing_list=sorted(contributing_list, key=lambda diagram: feynman_loop_number(diagram))
    return contributing_list

