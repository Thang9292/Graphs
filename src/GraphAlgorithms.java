import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * Implementation of various different graph algorithms.
 *
 * @author Thang Huynh
 * @version 1.0
 *
 */
public class GraphAlgorithms {

    /**
     * Performs a breadth first search (bfs) on the input graph, starting at
     * the parameterized starting vertex.
     *
     * When exploring a vertex, explore in the order of neighbors returned by
     * the adjacency list.The graph should be unmodified after this method terminates.
     *
     * @param <T>   the generic typing of the data
     * @param start the vertex to begin the bfs on
     * @param graph the graph to search through
     * @return list of vertices in visited order
     * @throws IllegalArgumentException if any input is null, or if start
     *                                  doesn't exist in the graph
     */
    public static <T> List<Vertex<T>> bfs(Vertex<T> start, Graph<T> graph) {
        //The Exceptions
        if (start == null || graph == null) {
            throw new IllegalArgumentException("The graph entered or the vertex entered was null");
        } else if (!graph.getVertices().contains(start)) {
            throw new IllegalArgumentException("The graph does not contain the starting vertex");
        }

        //Initialize Visited Set
        Set<Vertex<T>> visitedSet = new HashSet<>();
        //Initialize Queue
        Queue<Vertex<T>> theQueue = new LinkedList<>();
        //Initializes the adjacent list map
        Map<Vertex<T>, List<VertexDistance<T>>> adjList = graph.getAdjList();
        //Initialize the list of vertices we visited using BFS
        List<Vertex<T>> theReturningList = new ArrayList<>();

        //Adds start to visited Set
        visitedSet.add(start);

        //Enqueues start
        theQueue.add(start);

        //While the Queue is not empty
        while (!theQueue.isEmpty()) {

            //Dequeues and sets the current vertex
            Vertex<T> vertex = theQueue.peek();
            theQueue.remove();

            //Adds the vertex to our list that we will return
            theReturningList.add(vertex);

            //For all vertices (w) adjacent to our current vertex
            List<VertexDistance<T>> theList = adjList.get(vertex);
            for (VertexDistance<T> vertexDistance : theList) {
                Vertex<T> w = vertexDistance.getVertex();

                //If W is not in the visitedSet
                if (!visitedSet.contains(w)) {
                    //Mark W into the visited set
                    visitedSet.add(w);
                    //Enqueue w
                    theQueue.add(w);
                }
            }

        }
        //Returns the List of the vertices that we visited
        return theReturningList;
    }

    /**
     * Performs a depth first search (dfs) on the input graph, starting at
     * the parameterized starting vertex.
     *
     * When exploring a vertex, explore in the order of neighbors returned by
     * the adjacency list.
     *
     * This method is implemented recursively. The graph should be unmodified
     * after this method terminates.
     *
     * @param <T>   the generic typing of the data
     * @param start the vertex to begin the dfs on
     * @param graph the graph to search through
     * @return list of vertices in visited order
     * @throws IllegalArgumentException if any input is null, or if start
     *                                  doesn't exist in the graph
     */
    public static <T> List<Vertex<T>> dfs(Vertex<T> start, Graph<T> graph) {
        //The Exceptions
        if (start == null || graph == null) {
            throw new IllegalArgumentException("The graph entered or the vertex entered was null");
        } else if (!graph.getVertices().contains(start)) {
            throw new IllegalArgumentException("The graph does not contain the starting vertex");
        }

        //Initialize Visited Set
        Set<Vertex<T>> visitedSet = new HashSet<>();
        //Initializes the adjacent list map
        Map<Vertex<T>, List<VertexDistance<T>>> adjList = graph.getAdjList();
        //Initialize the list of vertices we visited using DFS
        List<Vertex<T>> theReturningList = new ArrayList<>();

        //Start of the recurse call
        rDFS(start, visitedSet, adjList, theReturningList);

        //Returns the List of the vertices that we visited
        return theReturningList;
    }

    /**
     * A helper method for doing recursive depth first search
     *
     * @param currentVertex     the current vertex we are looking at
     * @param visitedSet        the set of vertex we have marked as visited
     * @param adjList           the list of VertexDistance
     * @param theReturningList  the list containing all the vertices in the order we visit them
     * @param <T>               the generic typing of the data
     */
    public static <T> void rDFS(Vertex<T> currentVertex,
                                Set<Vertex<T>> visitedSet,
                                Map<Vertex<T>, List<VertexDistance<T>>> adjList,
                                List<Vertex<T>> theReturningList) {

        //Marks the current vertex as visited
        visitedSet.add(currentVertex);
        //Adds the current vertex as to our returningList
        theReturningList.add(currentVertex);

        //For all vertices (w) adjacent to our current vertex
        List<VertexDistance<T>> theList = adjList.get(currentVertex);
        for (VertexDistance<T> vertexDistance : theList) {
            Vertex<T> w = vertexDistance.getVertex();
            //If W is not in the visitedSet
            if (!visitedSet.contains(w)) {
                rDFS(w, visitedSet, adjList, theReturningList);
            }
        }
    }

    /**
     * Finds the single-source shortest distance between the start vertex and
     * all vertices given a weighted graph. We assume non-negative edge
     * weights
     *
     * Return a map of the shortest distances such that the key of each entry
     * is a node in the graph and the value for the key is the shortest distance
     * to that node from start, or Integer.MAX_VALUE (representing
     * infinity) if no path exists.
     *
     * The version of Dijkstra's where we use two
     * termination conditions in conjunction.
     *
     * 1) Check if all of the vertices have been visited.
     * 2) Check if the PQ is empty yet.
     *
     * @param <T>   the generic typing of the data
     * @param start the vertex to begin the Dijkstra's on (source)
     * @param graph the graph we are applying Dijkstra's to
     * @return a map of the shortest distances from start to every
     * other node in the graph
     * @throws IllegalArgumentException if any input is null, or if start
     *                                  doesn't exist in the graph.
     */
    public static <T> Map<Vertex<T>, Integer> dijkstras(Vertex<T> start,
                                                        Graph<T> graph) {
        //The Exceptions
        if (start == null || graph == null) {
            throw new IllegalArgumentException("The graph entered or the vertex entered was null");
        } else if (!graph.getVertices().contains(start)) {
            throw new IllegalArgumentException("The graph does not contain the starting vertex");
        }

        //Initialize the visitedSet
        Set<Vertex<T>> visitedSet = new HashSet<>();
        //Initialize the distance map
        Map<Vertex<T>, Integer> distanceMap = new HashMap<>();
        //Initialize the priority queue of vertex distance
        PriorityQueue<VertexDistance<T>> priorityQueue = new PriorityQueue<>();
        //Initialize the adjacent list map
        Map<Vertex<T>, List<VertexDistance<T>>> adjList = graph.getAdjList();

        //Sets the distance for start in the map to be 0
        distanceMap.put(start, 0);
        //For all vertices (v) in the graph (g), initialize distance of v to INF
        for (Vertex<T> vertex: graph.getVertices()) {
            //In this case, we are setting INF as the maximum value an int can hold
            if (!vertex.equals(start)) {
                distanceMap.put(vertex, 2147483647);
            }
        }

        //Enqueue (s, 0) into priority queue
        priorityQueue.add(new VertexDistance<T>(start, 0));

        //while PQ is not empty and VS is not full
        while (!priorityQueue.isEmpty() && visitedSet.size() < adjList.size()) {

            //(u, d) <-- PQ.dequeue() where u is the current vertex and d is the distance
            VertexDistance<T> vertexDistance = priorityQueue.peek();
            priorityQueue.remove();
            Vertex<T> currentVertex = vertexDistance.getVertex();
            int d = vertexDistance.getDistance();

            //if u is not visited in VS
            if (!visitedSet.contains(currentVertex)) {
                //Mark U as visited
                visitedSet.add(currentVertex);

                //Update the DistanceMap for u with new shortest path d
                //Only does this if the new distance is less than the old distance
                if (d < distanceMap.get(currentVertex)) {
                    distanceMap.put(currentVertex, d);
                }

                //for all (w, d2) adjacent to u and not visited in the visited set
                List<VertexDistance<T>> theList = adjList.get(currentVertex);
                for (VertexDistance<T> vertexDistance1 : theList) {
                    Vertex<T> w = vertexDistance1.getVertex();
                    int d2 = vertexDistance1.getDistance();

                    //W is not in the visited set
                    if (!visitedSet.contains(w)) {

                        //PQ.enqueue((w, d+d2))
                        priorityQueue.add(new VertexDistance<>(w, d + d2));
                    }
                }
            }
        }

        //Returns the distance map that contains the shortest path
        return distanceMap;

    }

    /**
     * Runs Kruskal's algorithm on the given graph and returns the Minimal
     * Spanning Tree (MST) in the form of a set of Edges. If the graph is
     * disconnected and therefore no valid MST exists, return null.
     *
     * Assume that the passed in graph is undirected. In this framework,
     * this means that if (u, v, 3) is in the graph, then the opposite edge
     * (v, u, 3) will also be in the graph, though as a separate Edge object.
     *
     * The returned set of edges should form an undirected graph. This means
     * that every time we add an edge to our return set, we should add the
     * reverse edge to the set as well. This is for testing purposes. This
     * reverse edge does not need to be the one from the graph itself; we can
     * just make a new edge object representing the reverse edge.
     *
     * Assume that there will only be one valid MST that can be formed.
     *
     * Kruskal's requires using a Disjoint Set which will keep track of which vertices are
     * connected given the edges in our current MST, allowing us to easily
     * figure out whether adding an edge will create a cycle. Refer
     * to the DisjointSet and DisjointSetNode classes for more information.
     *
     * Should NOT allow self-loops or parallel edges into the MST.
     *
     * By using the Disjoint Set provided, we can avoid adding self-loops and
     * parallel edges into the MST.
     *
     * @param <T>   the generic typing of the data
     * @param graph the graph we are applying Kruskals to
     * @return the MST of the graph or null if there is no valid MST
     * @throws IllegalArgumentException if any input is null
     */
    public static <T> Set<Edge<T>> kruskals(Graph<T> graph) {
        //The Exceptions
        if (graph == null) {
            throw new IllegalArgumentException("The graph entered was null");
        }
        //The graph is disconnected
        if (graph.getEdges() == null || graph.getEdges().size() == 0) {
            return null;
        }

        //Initialize DisjointSet, DS, with all vertices in G
        DisjointSet<Vertex<T>> disjointSet = new DisjointSet<>();

        //initialize MST EdgeSet, MST
        Set<Edge<T>> minSpanTree = new HashSet<>();

        //initialize PriorityQueue, PQ, with all edges in G using build heap
        PriorityQueue<Edge<T>> priorityQueue = new PriorityQueue<>(graph.getEdges());

        //while PQ is not empty and MST has fewer than n-1 edges
        while (!priorityQueue.isEmpty() && minSpanTree.size() < graph.getEdges().size() - 1) {

            //edge(u, v) <-- PQ.dequeue() where U and V are the vertices connected by that edge
            Edge<T> edge = priorityQueue.remove();
            Vertex<T> u = edge.getU();
            Vertex<T> v = edge.getV();

            //if u and v are not in the same cluster
            if (disjointSet.find(u) != disjointSet.find(v)) {
                //add edge(u, v) to MST
                minSpanTree.add(edge);

                //For Testing purpose only (adding in reverse edge)
                Edge<T> reverseEdge = new Edge<>(v, u, edge.getWeight());
                minSpanTree.add(reverseEdge);

                //merge u’s cluster with v’s cluster
                disjointSet.union(u, v);
            }
        }

        //Return the minimal spanning tree
        return minSpanTree;
    }
}
