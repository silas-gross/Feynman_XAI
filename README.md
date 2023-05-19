# Feynman_XAI
Final project for CUNY Graduate AI class. 
Using Feynman diagrams to set a renormalization scheme on a dataset and comparing the structure induced by the renomalization scheme to that of a bayes net learned from the data to improve infernce speeds and towards an exact mapping between renormalization and AI models in a similar way to DNNs as was shown by Mehta and Schwab (arXiv:1410.3831)

#Goals for this project 
##Overall
Broadly this project aims to solve a number of problems simultanteously. The First use of this is yet another Feynman diagram generator, a method to search for contributions to Scattering amplitudes, Willsonian and mass renomalization/running coupling renormalization schemes and implementing them in real time in the expansion, and establishing an equivelence with Bayes Net models.

##Short Term goals
For the purposes of the project deadline of 18 May 2022, the scope of this project is going to be limited to the first chunk of this project. 
The Goals are as follows 
1. Create a model that allows for numerical representation of Feynman diagrams with UV cut-off given from outside, assuming all exterior particles are at rest
2. Create an expansion method that allows for new vertices to be added to a Feynman diagram and recombination of the resultant particles onto other free particles.
3. Take in a Lagrangian from a GUI 
4. Read out vertices from Lagrangian
5. Calculate propagators from Lagrangian
6. Determine heuristic to expand vertices
7. Graphical representation of expansion as output
8. Build A\* search tree (potentially implement with LAO\* however the test Lagrangian is unlikely to need this)
9. Implement A\* with CSP, CSP constrains the expansions to fit within the RG Flow equation/ CS equation 
11. Calculate a stop point based on change in UV cutoff or running coupling change (finding fixed point)
12. Output scattering amplitudes for both methods

##Long term goals
Ideally the next step is to show that the running coupling approach yeilds the same structure in the Bayes net as the Willsonian.
Then if this can be established, aim to benchmark learning of this Feynman Diagram approach verus learning a Bayes net to see if either has any particular advantage, thus allowing for faster caluclation of one or the other.
In the feynman diagram represention, next steps could be making sure it accounts for all ways to recombine the particles, more general loops, allow for massive particles and non-scalars. 

#Aproach
This is a broad descriptive section, for a more precise version of the approach, please see the included control flow diagram.

The idea is that the first step is to expand out all vertex diagrams as the base level call of the Feynman Diagram method. Then, expansion is always limite to those diagrams that add loop level corrections to this diagram, has to have the same broad inflow and outflow particles to the diagram to be able to properly compare scattering amplitudes. 
Then the loop 0 search will select the largest change to UV correction from the diagrams and expands further along that line of interactions. Every step on the search then establishes new values of the UV cut-off and applies up the chain to establish values on further expansion, where expansion of a diagram is taking one additional vertex along a propagator line to generate further graphs. 

This approach is broadly just a graph expansion algorithm that has two characteristics assigned to each graph, and the second characteristic is the one being regulated from the UV divergence, but the first characteristic (the loop order) is teh characteristic that defines the depth of the search and the algorithm itself is limited by the CSP informed by UV cutoff. 

The goal state is where is defined as the state wherein a fixed point on the RG flow equation/ Callan-Syzemic equation (see comments in renormalization.py for exact formulation) is reached, or the requirements of the CSP have reached a state wherein it forces back-flow on the RG. 

So the flow of this is 

Get Largangian $\rightarrow$ Get Vertex diagrams $\rightarrow$ intiallized search with start value of UV cut off $\rightarrow$ establish channel of intrest $\rightarrow$ perform search on this, expansions having contributions limited by CSPs and informed by search $\rightarrow$ find solution, and compare the two UV cutoffs

#Challenges
Turns out doing obnoxiously hard math is even harder in the abstracting needed to represent this in the language of A\* search. Also divergent integrals turn out to be technically impossible to calculate and I am having to find some tricks to get around that. 
The Group representation that is needed is more difficult than anticpated, hence the narrowing of the project to its current goals, which is actually a limitation of procedural programing, as most packages in the broad space of Feynman Diagram representation and calulation are done in functional paradigms.

Also it is not clear from my reading the exact correlation between energy model Bayes nets and Bayes nets how we have defined them in the class, thus will need to be expanded on at later use.

#Usage
Usage is fairly simple, run python3 main.py and it will pop open a window to enter the lagrangian, and then proceed to run the search, outputing the found UV cutoff and scattering amplitudes when finished running, as well as a feynman diagram of the inital vertex and the final diagram  

#Data structure
So the main data structure for this project is an tuple that describes the state of a feynman diagram, including the number of particles, number of verticies, number of outgoing particles, loop order, scattering amplitude and a graph description of the particles. 
So, in form with () denoting tuples, [] denoting lists, {} denoting dictionaries
 	Total diagram:
	( [string form of incoming to outgoing particles], {all connections in the diagram}, scattering amplitude numerical value, number of incoming-outgoing particles, loop order of the diagram)

	Sub component 1: String form
	[ "Particles in $\rightarrow$ Particles out", [Particles in (as strings)], [Particles out (as strings)] ]

	Subcomponent 2: all connections in diagram 
	{node:(vertex contribution of this node, {connected node: connector})}

Subcomponent 2 gives a dictionary definition of the diagram, which can then be be visualized through the DrawDiagram() method in the Feynman_generator class. Ideally this method would be called as the diagrams are being parsed, giving a live status of the analysis, however doing such has been disabled in the current version, as I was not able to get the live update to work without blocking further running. Additionally, the current method does not account for direction of the graph which is an issue discussed below, but comes down to the way that the diagram dictionary is built, as the dictionary is technically doubling connections, connecting each node twice, i.e if the key node A gives a connection to node B, then key node B will give a connection to node A. This can be fixed when calling the AssignValues() method, however, as this method is a CSP, it does tend to add quite a bit of time to the computation. 

#Running Flow

For a detailed control flow diagram, see document "Control Flow.pdf" in this repository. The basic model is as follows 

main.py serves as the entry point to the program. This is where the GUI lives. main invokes an instance of the GUI, allowing the user to input the lagrangian over which we will run, as well as initial values for the masses of the particles associated with each field and the inital values of the coupling constants for each interaction term in the lagrangian. 

Then, a Feynman\_Search object in initialized with the Lagrangian as input, which in turn calls a Feynman\_Generator object in order to build diagrams. The Feynman\_Generator object by default initalizes with a set of base diagrams corresponing to each vertex given in the Lagrangian. That set of diagrams is passed first to the Feynman\_Search object as the base\_diagrams property, which is then passed back to main, and are passed to a list of "visual diagrams" which accepts and labels the diagrams to be displayed. 

Then, based on the Lagrangian and the base diagrams, main then calls the DecidePointFunction() method of the Feynman_Search object. This parses through the base diagrams, and searches the immediate children, i.e. tree diagrams with two vertices if possible or 1 loop corrections, of each base diagram to identify which one will generate the mazimum correction to each parameter immediately averaged over the number of children diagrams. This approach gives a decent proxy for benefit weighted by computational cost, approximating the calculations to fist order while accounting for the size of the tree that would need to be built in order to reach the same final goal, as the total final corrections on any one path will alway evaluate the same, but will require greater or fewer levels depending on the diagram-to-diagram corrections. Hence, if the correction is of order 1/o, where o is the original value of the coupling constant, from the first, fewer diagrams will be required to get to the asymptotic expansion (see report for explination of mathematics of asymptotic expansions and the limits of renomalization). 

Once this method has alight on one of the initial to final state diagrams from the base, that becomes the base for the construction of the A\* search algorithim. 
## AI algorithms

The main functionality of this project is based around 3 interacting AI systems, that all update eachother while running. The Feynman search method contains an A\* search that parses a tree formed by expandind the diagrams along vertex insertions, perfomed by the Feynman\_Generator object. This in turn, calls the ImposeConstraints() method. This method has two seperate running modes, Synthetic Coupling or Willsonian, corresponding to the two seperate renomalization approaches used in this project (see report for detailed description). If running in Synthetic Coupling mode, the system then performs a CSP search to impose a new value onto the coupling constant and the cutoff value for the integration with the constraint being informed by the current node being expanded by A\*, which in turn updates the heuristic definition within the A\* queue, by pushing the new values (and the requisite change in mass) onto the Feynman\_Generator object. In the Willisonian mode, the system simply imposes a function known as the beta function, which updates the coupling constant based on a fixed cutoff value parameter, thus avoiding the CSP, but missing out on the ability to properly coarse grain. 

The Feynman\_Generator object performs the actuall work of creating the child nodews and heuristic values, as well as creating the evaluation condition. The main AI system of this class is in the scattering amplitude calculation. Within the definition of how to calculate a scattering amplitude there is an ambiguity that arises as to how to assign momentum values to the internal lines of the diagram, and how to determine which particles are flowing into a vertex verus flowing out. Usually, when calculating by hand, one chooses this on a diagram by diagram basis given constraints that the momentum on each vertex must all add to 0, however, for this project, it becomes necessary to treat all momentums as free values and take a CSP approach by introducing the above constraint. Then, this is passed to the scattering amplitude integral, which depends on the coupling constants, masses and cutoff values, and produces the number which gives evaluation of the diagram to the A\* method that then uses this to update the ImposeConstraints() method for both running modes. 

Within the ImposeConstraints() method, the stop condition is encoded. This simply checks to see if the change in value has reached the point in the asymptotic expansion where the divergence from the proper value grows again, or if the correction is non-physical. Beyond that, the A\* search contians a further stop condition where in, if the change in scattering amplitude is lower than the hueristic for the next calculation even for a properly physical change, it is close enough to the correct answer that further searching would be inappropriate. 


