#This file will holdmath representation python the class to generate feynamn diagrams for searching 
#The generation will need to be a tuple to be hashable
#The hash then should be put in a map with a cost function to improve A* search
from math import comb
class FeynmanGenerator:
    def __init__(self, lagrangian, looporder=1):
        self.order=looporder
        self.vertices=self.GetVertices(lagrangian)
        self.propagator=self.GetPropagators(lagrangian)
        self.diagrams=self.GenerateDiagrams()
        self.scattering_amp=self.SumDiagrams()
        self.vertex_count=self.CountVertecies()
        self.diagram_output=self.GenerateOutput()
    def GenerateDiagrams():
        #this will generate diagrams up to loop order given 
        #if loop order=-1, generate until contribution to scattering amplitude is diminneshed by a power of lambda 
        #datastructure for a diagram needs to be 
        #{(string form of incoming to outgoing particles, incoming particles, outgoing particles), set of vertexes, (scattering amplitude numerical, scattering amplitude representational), ingoing-outgoing particles, loop order}
        #this will generate multiple diagrams filling the same set of charracteristics, thus, hash and every incidcdnce of a hash double used, add to the symmetry factor
        scattering_amplitude=0
        diagrams=list()
        hash_diagrams=list()
        vs=self.vertices
        ps=self.propagator
        for v in vs.keys(): #this has all of the single vertex configuration diagrams
            parts=v.split(",");
            unique_parts=dict()
            #need to loop over right and left hand symmetries
            #for all unique particles, there is nCk configurations for each 
            #config with k being the number of particles on the left hand 
            #the general for m many unique particles is a symmetry factor of 
            # nCk - (n-m)*min{k, (n-k)} 
            #so each set has an in-out of 
            #2k-n
            for i in parts:
                if i in unique_parts.keys():
                    unique_parts[i]+=1
                else:
                    unique_parts[i]=1
        #this has set up the number of particles at play and the number of unique particles
        #so now need to generate diagrams with the total number of particles and switch unquies between left and right
            n=len(parts)
            m=len(unique_parts.keys())
            for k in range(n):
                if k==0 or k==n:
                    continue
                sf=(n-m)*min(k, n-k)
                ndiags=comb(n,k) - sf
                for i in range(ndiags):
                    #need to select the particles incoming versus outgoing
                    l=list() #gives list of particles on the left
                    for p in range(k+1):
                        ix=p+i
                        ix=ix%n
                        l.append(parts[ix])
                    rparts=parts.copy()
                    for p in l:
                        rparts.remove(p)
                    sform=''.join(str(x) for x in l)
                    sform+="->"
                    sform=sform.join(str(x) for x in rparts)
                    d=((sform, l, rparts), {'1':[v, [{}]]}, (vs[v]/sf, "-i "+str(vs[v])+" 1/"+str(sf)), 2k-n, 0)
                    if not hash(d) in hash_diagrams:
                        diagrams.append(d)
                        hash_diagrams.append(d)
        #now I need to generate diagrams up to loop order l 
        # this is calculated with L=I-V+1
        # where I is the number of internal lines and V is number of vertices
        for d in diagrams:
            if d[0] in scattering_amplitude:
                scattering_amplitude[d[0][0]]+=d[2][0]
            else:
                scattering_amplitude[d[0][0]]+=d[2][0]
        #first we need to set the current value of the scattering amplitude of each in state versus outstate
        #now expand out each diagram from single vertexes to build a new set of diagrams
        all_added=False
        while all_added==False:
            one_branch=False
            for d in diagrams:
                d_add=self.GenerateNextOrder(d)
                for j in d_add:
                    if hash(j) in hash_diagrams:
                        d_add.remove(j)
                        continue
                    else:
                        hash_diagrams.append(hash(j))
                        one_branch=True
                        diagrams.append(j)
                        continue
            if !one_branch:
                all_added=True
                break
        if all_added==True:
            return diagrams
    def Generate_Next_order(self, diagram):
        d_add=self.ExpandDiagram(diagram)
        vs=self.vertices
        for i in d_add:
            if i[0][0] != d[0][0]:
                d_add.remove(i)
                second_order=self.ExpandDiagram(i) #allows to expand one extra vertex on the diagram to try to get to recombine to the desired state
                #this acutally needs to be refined to allow for recombination and the expansion will always just branch more out
                #can try to fix in expand diagrams
                for j in second_order:
                    if j[0][0] != d[0][0]:
                        continue
                    else:
                        d_add.append[j]
            else:
                if i[-1] > self.order:
                    d_add.remove(i)
                    continue
                sa=i[2][0]
                if self.order !=0:
                    cutoff= scattering_amplitude[i[0][0]]*pow(min(vs.values()), self.order) #limits contributions to those supressed by less than the loop order
                    if sa <= cutoff:
                        d_add.remove(i)
                        continue
                else:
                    second_order=self.ExpandDiagram(i)
                    ssa=0
                    for j in second_order:
                        ssa+=j[2][0]
                    first_ord_correction= sa/scattering_ampiltude[i[0][0]]
                    second_ord_correction=ssa/sa
                    if second_ord_correction <= first_ord_correction*min(vs.values()):
                        #this gives a text to see if the corrections are addign diminishing corrections
                        d_add.remove(i)
                        continue
        return d_add
    def ExpandDiagram(self, diagram):
        #this will take in the diagrams and apply the relevant propagators and vertices to get to a further stage
        vs=self.vertices
        ps=self.propagators
        dvs=diagram[1].copy()
        dvsk=dvs.keys()
        new_diagrams=list()
        dvsd=dvs.copy()
        for k in dvsk:
            vp=dvsk[k] #gives vertex and which it connects to 
            v=vp[0] #vertex
            pc=vp[1] #dictionary of connections and propagators that lead to them
            particles=v.split[","]
            for part in particles:
                if not part in pc.values():
                    #this section is specifically adding a particle that does not already appear in the list of particles
                    #will need to do the same for ones in, just broke it out here to be more clear
                    A=len(dvsk)
                    A+=1 #this would be to add a new node
                    pc[A]=part
                    dvsd[k]=[v,pc]
                    #This line takes the propagator in question and changes it such that it only goes to this one new propagator
                    for vn in vs:
                        if part in vn:
                            dvsd[A]=[vn, {k, part}]
                            passed_val=False
                            i=0
                            for oparts in vn:

                                if oparts==parts and not passed_val:
                                    passed_val=True
                                    continue
                                else:
                                    out_node="exterior_"+str(A)+"_"+str(i)
                                    i+=1
                                    dvsd[out_node]=[{}, {k, opart}]
                                    #not quite sure about the k here but this is the general idea
                                
                        
                            #this has added a new vertex to the graph
                            #this vertex is an exterior vertex
                else:
                    #need to allow for recombination on exterior lines
                    #for this I need very specific notation to exterior nodes
                    #have to add another type of key to my datastructure
        #recombine after expansion is done
        #now need to get propagators corresponding to the connecting arcs

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
        
        
