#This file will holdmath representation python the class to generate feynamn diagrams for searching 
#The generation will need to be a tuple to be hashable
#The hash then should be put in a map with a cost function to improve A* search
from math import comb
from scipy import integrate
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
                    out=dict() #this is the dictionary of the connections to the exterior
                    for j in range(len(parts)):
                        out_node="exterior_"+str(j)
                        out[out_node]=parts[j]
                    d=((sform, l, rparts), {'1':[v, out]}, (vs[v]/sf, "-i "+str(vs[v])+" 1/"+str(sf)), 2k-n, 0)
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
    def ExpandDiagram(self, diagram, leg):
        #this will take in a diagram and expand the vertex along a specific propagator
        #redoing the heuristic so that it just accounts for the mass on the propagator as a proxy for a scaling on the SA
        #so it will be average of coupling constant over the vertexes including the propagator /m^2
        #that is a good proxy to the scattering amplitude with out calcuating
        #so this is what the "expand children" will be 
        vs=self.vertices
        ps=self.propagators

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
        #each set will be a vertex and the propagators flowwing out
        return sa
    def heuristic(self, diagram, propagator):
        #the heuristic here is given by h=1/SA*(average lambda)/m_prop^2+1-(in+out)/total lines 
        #this last part is a normalized proxy for the loop order
        #all this is saying is really just contributions to scattering amplitude and order of divergence
        avg_constant=0
        total_vertexs=0
        vertexs=self.vertices.copy()
        in_out=len(diagram[0][1])+len(diagram[0][2]) #this is base level
        n_lines=sum([len(diagram[1][x]) for x in diagram[1].keys])
        in_out_new_contribution=0
        for vk in vertexs.keys():
            if propagator in vk:
                avg_constant+=vertexs[vk]
                total_vertexs+=1
                in_out_new_contribution+=len(vk.split(","))-1
        if total_vertexs !=0:
            avg_constant=float(avg_constant)/total_vertexs
            in_out_new_contribution=in_out_new_contribution/total_vertexs
        SA=diagram[2][0] #current scattering amplitude
        m=self.propagators[p][0] #mass of the propagator
        starting_coupling_constant=self.CountVertices(diagram)[1]
        #so I need count vertices to return a list of the form 
        #(count of vertexs, product of coupling constants, list of vertexs)
        h=1/SA*(starting_coupling_constant*avg_constant)/pow(m,2)
        n_lines+=in_out_new_contribution
        in_out+=in_out_new_contribution
        h+= 1-in_out/n_lines
        return h
        

             
        
        
