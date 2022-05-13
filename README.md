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
Ideally the nextstep is to show that the running coupling approach yeilds the same structure in the Bayes net as the Willsonian.
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

