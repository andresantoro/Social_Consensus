FastConsensusSolver
===================

Class to solve a continuous opinion formation problem. The project is built with [pybind11](https://github.com/pybind/pybind11).

Installation
------------

 `pip install ./FastConsensusSolver`

or developer mode :

 `pip install -e ./FastConsensusSolver`

Example
-------

### Consensus on quenched network

In this example, we define the structure as a complete graph and we assign to every node the same influence. We use the solver to reach a consensus on the opinion for each node, then we draw the temporal evolution of the process.

```python
from FastConsensusSolver import QuenchedSolver
import matplotlib.pyplot as plt
import numpy as np

#define the system
network_dict = {0: [1,2,3,4], 1: [0,2,3,4], 2: [1,0,3,4], 3: [1,2,0,4], 4: [1,2,3,0]}
influence_dict = {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5}
eta = 0.5
seed = 42

#define the solver with initial random state
cs = QuenchedSolver(network_dict, influence_dict, eta, seed)

#reach consensus
tol = 0.005
cs.reach_consensus(tol)

#get the raw data for temporal evolution of state
initial_state = cs.get_initial_state_vector()
history_vector = cs.get_history_vector()

#get the variation from the history and update state
node_vector = np.arange(0,4)
time_vector = np.arange(len(history_vector)+1)
temporal_state_dict = {node: [opinion] for node, opinion in zip(node_vector, initial_state)}
for t in time_vector[:-1]:
    for node, temporal_opinion_vector in temporal_state_dict.items():
        if history_vector[t][0] == node:
            temporal_opinion_vector.append(temporal_opinion_vector[-1]+history_vector[t][1])
        else:
            temporal_opinion_vector.append(temporal_opinion_vector[-1])

#draw the opinion of each node
for node in node_vector:
    plt.plot(time_vector, temporal_state_dict[node])
plt.xlabel(r"Time")
plt.ylabel(r"Opinion")
plt.show()
```

### Consensus on annealed network

In this example, there is no structure. We assign to every node a different influence and a different priority to speak. We use the solver to reach a consensus on the opinion for each node, then we draw the temporal evolution of the process.

```python
from FastConsensusSolver import AnnealedSolver
import matplotlib.pyplot as plt
import numpy as np

#define the system
N = 20
priority_dict = {i : i for i in range(N)}
influence_dict = {i: (i+1)/N for i in range(N)}
eta = 0.5
seed = 42

#define the solver with initial random state
cs = AnnealedSolver(priority_dict, influence_dict, eta, seed)

#reach consensus
tol = 0.005
cs.reach_consensus(tol)

#get the raw data for temporal evolution of state
initial_state = cs.get_initial_state_vector()
history_vector = cs.get_history_vector()

#get the variation from the history and update state
node_vector = np.arange(0,N)
time_vector = np.arange(len(history_vector)+1)
temporal_state_dict = {node: [opinion] for node, opinion in zip(node_vector, initial_state)}
for t in time_vector[:-1]:
    for node, temporal_opinion_vector in temporal_state_dict.items():
        if history_vector[t][0] == node:
            temporal_opinion_vector.append(temporal_opinion_vector[-1]+history_vector[t][1])
        else:
            temporal_opinion_vector.append(temporal_opinion_vector[-1])

#draw the opinion of each node
for node in node_vector:
    plt.plot(time_vector, temporal_state_dict[node])
plt.xlabel(r"Time")
plt.ylabel(r"Opinion")
plt.show()
```

