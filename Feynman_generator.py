#This file will hold the class to generate feynamn diagrams for searching 
#The generation will need to be a tuple to be hashable
#The hash then should be put in a map with a cost function to improve A* search
class FeynmanGenerator:
    def __init__(self, lagrangian, looporder=1):
        self.order=order
        self.vertices=self.GetVertices(lagrangian)
        self.propagator=self.GetPropagators(lagrangian)
        self.diagrams=self.GenerateDiagrams(lagrangian)
        self.scattering_amp=self.SumDiagrams()
        self.vertex_count=self.CountVertecies()
        self.diagram_output=self.GenerateOutput()
    def GenerateDiagrams(lagrangian):
        #this will generate diagrams up to loop order given 
        #if loop order=-1, generate until contribution to scattering amplitude is diminneshed by a power of lambda 
        #datastructure for a diagram needs to be 
        #{set of vertexes, scattering amplitude, ingoing-outgoing particles, loop order}
        #this will gnerate multiple diagrams filling the same set of charracteristics, thus, hash and every incidcdnce of a hash double used, add to the symmetry factor
        return diagrams
    def GenerateOutput(self):
        #want to output a list of diagrams with associated costs
        #do this as a dictionary on the hashed tuple
        do=dict()
        a=self.scattering_amp
        o=self.order
        if o== -1:
            o=1
        for d in self.diagrams:
            hv=hash(d)
            c=(1-d[1]/a)+ d[2] + 1/o*(len(d[0])+d[3]) #this is a normalized-ish cost function
            do[hv]=c
        return do
    def SumDiagrams(self):
        #This will sum up all the present diagrams in the list of diagrams
        diag=self.diagrams
        sum_of_amplitudes=0
        for l in diag:
            sum_of_amplitudes+=l[1]
        return sum_of_amplitudes
    def CountVertices(self):
        diag=self.diagrams
        vertexs=dict()
        for d in diag:
            vl=d[1]
            for v in vl:
                if v in vertexs:
                    vertexs[v]+=1
                else:
                    vertexs[v]=1
        return vertexs
    def GetVertices(self, lagrangian):
        #read in the interaction parts of the lagrangian and give a dictionary of vertecies with correponding coupling constants
        #maybe hash down the road to compare models?
        vs=dict()
        return vs
    def GetPropagators(self, lagrangian)
        #calculate the propagator for each particle that we need to account for
        #return dictionary 
        ps=dict()
        return ps
    def CalculateScatteringAmplitude(self, diagram):
        #this method takes in a digram in the form of 
        #{ordered set of vertecies,connecting propagators, exterior label contirbutions}
        #each set will be a vertex and the propagators flowwing out
        sa=diagram[1]
        for vc in diagram[0]:
            #each element of diagram 
            sa+=self.verticies(vc[0])
            for p in vc[1]:
                sa+=self.propagators(p)
        return sa
        
        
