import networkx as nx
import numpy as np
import json
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.animation import FuncAnimation

#sys.argv[1] -> file containing the network structure
#sys.argv[2] -> file containing the distribution of opinion over time


#Function that loads the graph structure from the file "filename" into a dictionary
#with the format: [key:nodes,values:first_neighbours]
def load_network_fromfile(filename):
    json_data=open(filename).read()
    graph=json.loads(json_data)
    graph={int(k):[int(i) for i in v] for k,v in graph.items()}
    return(graph)

#Function that creates the networkx structure from the edge list dictionary. It returns the graph G as output
def create_graph_fromedgelist(graph_structure):
    G=nx.Graph()
    for key in graph_structure.keys():
        for l in graph_structure[key]:
            G.add_edge(key,l)
    return(G)


# generate graph
graph_structure= load_network_fromfile(str(sys.argv[1]))
G= create_graph_fromedgelist(graph_structure)
color_map=np.loadtxt(str(sys.argv[2]))

color=[]
for i,line in enumerate(color_map):
	if i % 100 == 0:
		color.append(line)



# draw the topology of the graph, what changes during animation
# is just the color
pos = nx.spring_layout(G,k=0.25,iterations=20)
nodes = nx.draw_networkx_nodes(G,pos)
edges = nx.draw_networkx_edges(G,pos)

# nx.draw_networkx_nodes(G, pos, nodelist=nodes, \
#     node_color='blue', node_shape='o')
# nx.draw_networkx_nodes(G, pos, nodelist=nodes, \
#     node_color='purple', node_shape='s')
plt.axis('off')


# pass frames to funcanimation via update function
# this is where I get stuck, since I cannot break
# out of the loop, neither can I read every array of
# the ndarray without looping over it explicitly
def update(i):
    # for i in range(len(frame)):
    # instead of giving frame as input, if I randomly generate it, then it works
    nc = color[i]
    nodes.set_array(nc)
    return nodes,




# output animation; its important I save it
#
#fig = pl.figure(1)
fig = plt.gcf()
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = plt.colorbar(sm)



ani = FuncAnimation(fig, update, interval=20, frames=range(len(color)), blit=True)
ani.save('crap.gif', writer='imagemagick',  savefig_kwargs={'facecolor':'white'}, fps=5)





# graph_structure= load_network_fromfile(str(sys.argv[1]))
# G= create_graph_fromedgelist(graph_structure)


# pos=nx.spring_layout(G)
# color_map=np.loadtxt(str(sys.argv[2]))
# # for i in range(0,len(color_map)):
# # 	print(i)
# # 	nx.draw(G,node_color=color_map[i,:], vmin=0, vmax=1,pos=posix)
# # 	plt.show()
