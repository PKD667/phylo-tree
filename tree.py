from compare import taxon_distances
import numpy as np

from treelib import Tree
from typing import Optional, List, Dict, Any

import matplotlib.pyplot as plt



class PhyloTree:
    def __init__(self, name: str, parent: Optional['PhyloTree'] = None):
        """
        Initialize a PhyloTree node.
        
        :param name: Name of the taxon.
        :param parent: Parent node of the current node.
        """
        self.name = name
        self.parent = parent
        self.children: List['PhyloTree'] = []

    def get_child(self, taxon: str) -> Optional['PhyloTree']:
        """
        Get the child node corresponding to the given taxon.
        
        :param taxon: Name of the taxon to search for.
        :return: The child node if found, otherwise None.
        """
        for child in self.children:
            if child.name == taxon:
                return child
        return None

    def add_child(self, child: 'PhyloTree'):
        """
        Add a child tree to the current node.
        
        :param child: The child tree node to be added.
        """
        self.children.append(child)

    def visualize(self):
        """
        Visualize the phylogenetic tree using the treelib library.
        """
        tree = Tree()
        tree.create_node(self.name, self.name)  # Root node
    
        def add_children(node: 'PhyloTree', parent_name: str):
            for child in node.children:
                tree.create_node(child.name, child.name, parent=parent_name)
                add_children(child, child.name)
    
        add_children(self, self.name)
        print(tree.show(stdout=False))
        

def print_dists(taxons,distances):
    # print the taxons
    print("  ",end=" ")
    for i in range(0,len(taxons)):
        print(taxons[i],end="     ")
    print()

    for i in range(0,len(distances)):
        print(taxons[i],end=" ")
        for j in range(0,len(distances[i])):
            print(f"{distances[i][j]:.1f}",end=" ")
        print()

def graph_dists(taxons,distances):

    # create a heatmap of the distances
    fig, ax = plt.subplots()
    cax = ax.matshow(distances, cmap='hot')

    # Add colorbar
    fig.colorbar(cax)

    # plot the taxons
    ax.set_xticks(range(0,len(taxons)))
    ax.set_yticks(range(0,len(taxons)))

    ax.set_xticklabels(taxons)
    ax.set_yticklabels(taxons)

    plt.show()



# return for each index the order of the distances
def get_order_table(distances):
    order = []
    for i in range(0,len(distances)):
        order.append(np.argsort(distances[i]))

    # remove the first element of each array
    for i in range(0,len(order)):
        order[i] = order[i][1:]

    return order

# here we have the algorithm finding the root node of the tree
def get_parent(distances):
    dmin = np.inf
    dmin_i = None
    for i in range(0,len(distances)):
        if sum(distances[i]) < dmin :
            dmin = sum(distances[i])
            dmin_i = i
    return dmin_i

# algorithm to compute the children of a node
def get_children(order,parent):
    children = []

    # Loop through the order table
    for i in range(0,len(order)):

        # If the node is the parent, skip it
        if i == parent:
            continue
        
        # Loop through the order of the node
        for j in range(0,len(order[i])): 
            
            # If we found another children closer than the parent, 
            # It means that the node is deeper in the tree
            if order[i][j] in children:
                break
            # if we found the parent, we add the node to the children
            elif order[i][j] == parent:
                children.append(i)
                break
            
            # Here we are gonna check if the node is deeper in the tree
            # If the node preceding the parent in the order is deeper, 
            # it means that the preceding node is a deep node, and that our node
            # might be a child directly of the parent

            # K is the index of the node preceding the parent in the order
            k = order[i][j]
            is_deeper = False

            # Loop through the order of the node preceding the parent (k)
            for l in range(0,len(order[k])):

                # If we found the parent, it means that the node [ i ] is deeper in the tree
                if order[k][l] == parent:
                    # We break the loop to avoid adding the node [ i ] to the children
                    break
                
                # if we found our node[ i ] in the order of the preceding node[ k ]
                # it means that the node[ k ] is deeper in the tree
                if order[k][l] == i:
                    is_deeper = True
                    break

            # If the node [ k ] is deeper in the tree, we continue to the next node
            if not is_deeper:
                break

    return children

# get new distances
def get_new_distances(distances,taxons,children,parent):
    children_data = []
    for child in children:
        child_data = {
            "taxons" : [],
            "distances" : []
        }

        sub_indices = []

        # add all the node closer to child than to the other children
        for i in range(0,len(distances)):

            # if the node is the parent or a child, we skip it
            if i == parent or i in children:
                continue

            # Reference distance
            # distance from node[ child ] to node[ i ]
            dref = distances[i][child]
            is_closest = True

            for c in children:

                #if another child is closer to the node[ i ] than the node[ child ]
                # we skip the node[ i ]
                if distances[i][c] < dref:
                    is_closest = False
                    break
            
            # if the node[ i ] is closer to the node[ child ] than to the other children
            # we add the distance to the new distances
            if is_closest and distances[i][parent] > dref:
                sub_indices.append(i)
        
        # convert the sub indices we gathered to the new distances, and taxons.
        # in the new distances we only keep the distances between the child and the sub indices nodes

        for i in sub_indices:
            child_data["taxons"].append(taxons[i])
            child_data["distances"].append([distances[i][j] for j in (sub_indices + [child])])

        # then add the taxon and data from the child itself
        child_data["taxons"].append(taxons[child])
        child_data["distances"].append([distances[child][j] for j in (sub_indices + [child])])
        

        children_data.append(child_data)

    return children_data
        
    
# Build a tree using a parent specie.
def build_tree_dist(taxons,distances):

    parent = get_parent(distances)
    tree = PhyloTree(taxons[parent],parent=None)

    order = get_order_table(distances)

    children = get_children(order,parent)

    print(f"Parent: {taxons[parent]}")
    print(f"Children: {children}")

    if len(distances) == len(children) + 1:
        for c in children:
            ctree = PhyloTree(taxons[c])
            tree.add_child(ctree)
        
        return tree

    new_data = get_new_distances(distances,taxons,children,parent)

    for i in range(0,len(children)):
        d = new_data[i]["distances"]
        t = new_data[i]["taxons"]
        print_dists(t,d)

        child_tree = build_tree_dist(t,d)

        tree.add_child(child_tree)

    tree.visualize()

    return tree

def build_tree(taxons):
    distances = taxon_distances(taxons)
    return build_tree_dist(taxons,distances)
    





