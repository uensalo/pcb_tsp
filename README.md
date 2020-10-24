# pcb_tsp

This repository contains a project I did during my undergraduate studies which involved formulating the drill-head movement of a PCB drilling machine as an optimization problem. The problem is that we'd like to drill some predetermined positions on a PCB; however there are soldered electrical components on the board over which the drill head may not pass through. As we are efficient creatures, we'd like to save time and energy by minimizing the movement of the drill head while not going over any of the forbidden regions. We also assume that we are allowed to visit a hole to be drilled only once, to keep things similar to the travelling salesman problem.

The problem was solved by converting a discrete Cartesian grid into a graph, then running DFS to see if a node is reachable from one node to another, while also running Floyd-Warshall to see if the path through which this node is reachable is the path with the shortest Manhattan distance in MATLAB. Then, all the information was fed to a CPLEX program, which solves a TSP given the graph. I did the visualisation in Processing.

Mathematical derivations and relevant definitions are located under paper.pdf. MATLAB files are for preprocessing (grid to graph conversion and graph-weighing), and CPLEX project files are for solving the corresponding TSP.

This was a group project, and I'd like to thank them even for the smallest contribution they've made.
