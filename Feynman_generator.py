#This file will holdmath representation python the class to generate feynamn diagrams for searching 
#The generation will need to be a tuple to be hashable
#The hash then should be put in a map with a cost function to improve A* search
from math import comb
from copy import deepcopy
import sympy as smb
from random import choice
import networkx as nx
import matplotlib.pyplot as plt
import sys
import builtins

class FeynmanGenerator:
    def __init__(self, lagrangian, cutoff, looporder=1):
        #by default, evaluate the diagrams to loop order 1 
        #each time the generator is called, loop order increases by 1
        self.order=looporder
        self.cutoff=cutoff
        self.diagram_list=[]
        self.vertices=self.GetVertices(lagrangian)
        self.propagators=self.GetPropagators(lagrangian)
        self.diagrams=self.GenerateDiagrams()
        self.exterior=[]
#        self.diagram_list=list() #list of the hashs of diagrams
    def UpdateCutoffandVertex(self, new_cutoff, newccs):
        #Feynman_Search class updates the coupling constant and vertices 
        #this is the renormalization, which then changes the  integrals for the scattering amplitude
        #RG flow equations mean that the amplitude will be independent of renomalization
        #however, the cutoff value will potentially factor into the calculation
        #and the renormalization parameters are required in order to perform integrals
        self.cutoff=new_cutoff
        self.vertices=newccs
    def GenerateDiagrams(self):
        #this will generate diagrams up to loop order given 
        #if loop order=-1, generate until contribution to scattering amplitude is diminished by a power of the cutoff
        #callin loop order=-1 will cause a run time on order of a day 
        
        #datastructure for a diagram needs to be 
        #{(string form of incoming to outgoing particles, incoming particles, outgoing particles),
        #set of vertexes, scattering amplitude numerical,ingoing-outgoing particles, loop order}
        #this will generate multiple diagrams filling the same set of charracteristics,
        #thus, hash the diagram, as each component does not uniquely determine the diagram
        #rather it is the total diagram
        scattering_amplitude=0
        #scattering_amplitude is calculated later 
        #to calculate here would be too slow
        #rather, it is better to call as needed
        diagrams=list()
        hash_diagrams=self.diagram_list
        vs=self.vertices
        for v in vs.keys(): #this has all of the single vertex configuration diagrams
            parts=v.split(", ");
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
                if sf==0:
                    sf=1
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
                    sform=str()
                    for x in l: sform+=' '+str(x)
                    sform+="->"
                    for x in rparts: sform+=' '+ str(x)
                    if len(rparts) == 0: sform+=" vacuum"
                    out=dict() #this is the dictionary of the connections to the exterior
                    out_cons=dict()
                    for j in range(len(parts)):
                        out_node="out_"+str(j)
                        out[out_node]=parts[j]
                        out_cons[out_node]=[0, {'1':parts[j]}]
                        #connections between the out vertex and the node
                        #stucture of the dictionary node:(vertex contribution of this node, {connected node: connector})
                    out_cons['1']=[v, out]
                    r=2*k-n
                    d=[[sform, l, rparts], out_cons, vs[v]/sf,r , scattering_amplitude]
                    if not hash(str([*d])) in hash_diagrams:
                        diagrams.append(d)
                        hash_diagrams.append(hash(str([*d])))
        #now need to generate diagrams up to loop order L 
        # this is calculated with L=I-V+1
        # where I is the number of internal lines and V is number of vertices
        #first we need to set the current value of the scattering amplitude of each in state versus outstate
        #now expand out each diagram from single vertexes to build a new set of diagrams
        return diagrams
    def GenerateNextOrder(self, diagram_a):
        #need to figure out if I actually want to use this, if so, need to update to match with current Expand
        #what this is doing is actulalay generating the queue
        legs=list()
        prop_and_nodes=list()
        children=list()
        diagram=dict() 
        if type(diagram_a) is list:
            diagram=diagram_a[1]
        elif type(diagram_a) is dict:
            diagram=diagram_a
        for n in diagram.keys():
            for m in diagram[n][1].keys():
                leg=str(n)+":"+str(m)
                lgr=str(m)+":"+str(n)
                found=(leg in legs) or (lgr in legs)
                if not found:
                    legs.append(leg)
                    prop_and_nodes.append([diagram[n][1][m], n, m])
        for i in range(len(legs)):
            leg=legs[i]
            n=prop_and_nodes[i][1]
            m=prop_and_nodes[i][2]
    #        print(diagram)
            child_diagrams=self.ExpandDiagram(diagram, leg, n, m) 
            try:        #not sure why but this keeps giving an Nonetype object
     #           print(type(child_diagrams))
                if len(child_diagrams) == 0:
                    continue
                for c in child_diagrams:
                    children.append(self.GenerateOutput(c, prop_and_nodes[i]))
            except:
                children.append(self.GenerateOutput(diagram_a, prop_and_nodes[i]))
        return children
    def VT(self, particle_list):
        real_part=0
        complex_part=0
        for p in particle_list:
            if "phi" in str(p):
                real_part+=1
            if "psi" in str(p):
                complex_part+=1
        return [real_part, complex_part]
    def Same(self, particle_list_1, particle_list_2):
        none_diff=True
        for j in particle_list_1:
            if not j in particle_list_2:
                none_diff=False
        return none_diff

    def ExpandDiagram(self, diagram, leg, from_node, to_node, second=False, old_ext=[]):
        #this will take in a diagram and expand the vertex along a specific propagator
        #redoing the heuristic so that it just accounts for the mass on the propagator as a proxy for a scaling on the SA
        #so it will be average of coupling constant over the vertexes including the propagator /m^2
        #that is a good proxy to the scattering amplitude with out calcuating
        #so this is what the "expand children" will be 
        vs=self.vertices
        #just take in the graph form of the diagram
        out_diagrams=list()
        external_particles=self.exterior
        diagrams_hash=self.diagram_list
        n_out=len(external_particles)
        for v in vs:
            particles_in_vertex=v.split(", ")
            if leg in particles_in_vertex:
                new_external=list()
                remaining_particles=deepcopy(particles_in_vertex)
                remaining_particles.remove(leg)
                diagram_to_expand=deepcopy(diagram)
                connections=diagram_to_expand[from_node][1]
                #the idea is that you want to insert another node in between the two 
                if to_node in connections.keys(): 
                    try:
                        del connections[to_node]
                        del diagram_to_expand[to_node][1][from_node]
                    except:
                        print("from node " +str(from_node)+", to node "+str(to_node))
                else:
                    continue
                    #if the two nodes are not connected, the expansion does not make sense
                new_node=len(diagram_to_expand.keys())+1-n_out
                #create a new node that will be connected to the "from" node by the original particle
                connections[new_node]=leg
                diagram_to_expand[new_node]=[0, dict()]
                diagram_to_expand[new_node]=['',{from_node:leg}]
                #now we connect the new node to the old to node and allow for any new external particles
                #need to account for having already connected one of the particles in the vertex back
                #this is accomplished ysing the "remainging particles"
                for p in remaining_particles:
                    try:
                        diagram_to_expand[new_node][1][to_node]=p
                        diagram_to_expand[to_node][1][new_node]=p
                    except:
                        diagram_to_expand[new_node]=['', {from_node:leg, to_node:p}]
                        diagram_to_expand[to_node]=['', {new_node:p}]
                    print(diagram_to_expand[new_node])
                    new_particles=deepcopy(remaining_particles)
                    new_particles.remove(p)
                    for p1 in new_particles:
                        n_out+=1
                        out_node="out_"+str(n_out)
                        diagram_to_expand[new_node][1][out_node]=p1
                        diagram_to_expand[out_node]=['', {new_node:p1}]
                    #now we have connected all the nodes together and need to work on combining 
                    diagram_to_expand=self.ReadVertices(diagram_to_expand)
                    #this just sets the 0 component of the nodes properly for easy comparison/to get diagram working correctly
                    if diagram_to_expand[new_node][0] not in vs: 
                        print("bad vertex " +str(diagram_to_expand[new_node][0]))
                        particles=diagram_to_expand[new_node][0].split(' ')
                        if self.VT(particles) != self.VT(particles_in_vertex):
                            continue #makes sure new vertex is valid
                        else:
                            print("had weird ordering")
                            diagram_to_expand[new_node][0]=v
                    recombine=deepcopy(diagram_to_expand)
                    recombine=self.RecombineExteriors(recombine)
                    rc_new_external=list()
                    for x in diagram_to_expand.keys():
                        if "out_" in str(x):
                            if len(diagram_to_expand[x][1])==1:
                                new_external.append(diagram_to_expand[x][1].values())
                    print("new externals: "+str(new_external))
                    if len(recombine)==0: recombine=[deepcopy(diagram_to_expand)]
                    for x in recombine[0].keys():
                        if "out_" in str(x):
                            if len(recombine[0][x][1])==1:
                                rc_new_external.append(recombine[0][x][1].values())
                    if self.Same(rc_new_external, external_particles) or self.Same(new_external, external_particles): print("found a good val")
                    if not self.Same(new_external, external_particles):
                        if second==False:
                            for m in diagram_to_expand.keys():
                                for n in diagram_to_expand[m][1].keys():
                                    ds=self.ExpandDiagram(diagram_to_expand, diagram_to_expand[m][1][n], m, n, True, external_particles)
                                    for dss in ds:
                                        nep=list() 
                                        print("2nd order expanded " +str(dss))
                                        for node in dss.keys():
                                            if "out" in str(node):
                                                nep.append(list(dss[node][1].values())[0])
                                        if self.Same(nep,external_particles):
                                            tdss=str([*dss])
                                            if hash(tdss) in diagrams_hash:
                                                continue
                                            else:
                                                good=True
                                                for x in dss.keys():
                                                    if dss[x][0] not in vs:
                                                        good=False
                                                if good:
                                                    out_diagrams.append(dss)
                                                    diagrams_hash.append(tdss)
                                        else:
                                            continue
                    if not self.Same(rc_new_external, external_particles):
                        if second==False:
                            for m in recombine[0].keys():
                                for n in recombine[0][m][1].keys():
                                    ds=self.ExpandDiagram(recombine[0], recombine[0][m][1][n], m, n, True, external_particles)
                                    for dss in ds:
                                        nep=list()
                                        for node in dss.keys():
                                            if "out" in str(node):
                                                nep.append(list(dss[node][1].values())[0])
                                        if self.Same(nep,external_particles):
                                            tdss=str([*dss])
                                            if hash(tdss) in diagrams_hash:
                                                continue
                                            else:
                                                good=True
                                                for x in dss.keys():
                                                    if dss[x][0] not in vs:
                                                        good=False
                                                if good:
                                                    out_diagrams.append(dss)
                                                    diagrams_hash.append(tdss)
                                        else:
                                            continue
                    d=str(diagram_to_expand)

                    if hash(d) in diagrams_hash:
                        continue
                    else:
                        good=True
                        for x in diagram_to_expand.keys():
                            if diagram_to_expand[x][0] not in vs:
                                good=False
                        if good:
                            out_diagrams.append(diagram_to_expand)
                            diagrams_hash.append(d)
                    d=str(recombine)

                    if hash(d) in diagrams_hash:
                        continue
                    else:
                        good=True
                        for rc in recombine:
                            for x in rc.keys():
                                if rc[x][0] not in vs:
                                    good=False
                            if good:
                                out_diagrams.append(rc)
                                diagrams_hash.append(d)
                    #so now I've added all the diagrams 
                #now need to put in recombinations
            

            
                    to_recombine=True
                    dl=[diagram_to_expand]
                    while to_recombine:
                        dl=[self.RecombineExteriors(dx) for dx in dl]
                        dl=[x for y in dl for x in y]
                        if len(dl)==0:
                            to_recombine=False
                            break
                        else:
                            for dd in dl:
                                d=str(dd)
                                new_external=list()
                                for node in dd:
                                    if "out" in str(node):
                                        new_external.append(list(dd[node][1].values())[0])
                                if not self.Same(new_external, external_particles):
                                    pass
                                if hash(d) in diagrams_hash:
                                    continue
                                else:
                                    good=True
                                    for x in diagram_to_expand.keys():
                                        if diagram_to_expand[x][0] not in vs:
                                            good=False
                                    if good:
                                        out_diagrams.append(dd)
                                        diagrams_hash.append(d)
        if len(out_diagrams) == 0: out_diagrams=[diagram]
        print(len(out_diagrams))
        for odo in out_diagrams:
            if type(odo) is list:
                out_diagrams.remove(odo)
                for o in odo:
                    out_diagrams.append(o)
            else:
                continue
        #print(type(out_diagrams[0]))
        return out_diagrams


    def RecombineExteriors(self, diagram):
        #recombines exterior propagators
        #this keeps the number of incomming and outgoing particles stable
        #and thus allows us to look at the same process post selection
        #selection is done at vertex level, performed in Feynman_search
        k=[]
        for x in diagram.keys():
            if "out" in str(x):
                k.append(x)
        diagrams_out=list()
        dummy_diagram=deepcopy(diagram)
        d_hash=list()
        for i in range(len(k)):
            for j in range(len(k)):
                if i==j: 
                    continue
                else:
                    if not k[i] in dummy_diagram.keys() or not k[j] in dummy_diagram.keys():
                        continue
                    if  list(dummy_diagram[k[i]][1].values())[0]==list(dummy_diagram[k[j]][1].values())[0]:
                        c=list(dummy_diagram[k[j]][1].keys())[0]
                        b=list(dummy_diagram[k[i]][1].keys())[0]
                        p=list(dummy_diagram[k[i]][1].values())[0]
                        dummy_diagram[b][1][c]=p
                        dummy_diagram[c][1][b]=p
                        del dummy_diagram[k[j]]
                        del dummy_diagram[k[i]]
                        
                        d=str(dummy_diagram)
                        if hash(d) in d_hash:
                            continue
                        else:
                            diagrams_out.append(dummy_diagram)
                            d_hash.append(d)
        return diagrams_out

    def ReadVertices(self, diagram):
        for v in diagram:
            #print(diagram[v])
            legs=diagram[v][1].values()
            p=0
            for l in legs:
                if p == 0: 
                    p=''
                else:
                    p+=' '
                p+=str(l)
            diagram[v][0]=p
        return diagram
    def GenerateOutput(self, diagram, prop_and_nodes):
        #need to figure out what this should be
        output=[diagram, prop_and_nodes, self.heuristic(diagram, prop_and_nodes[0])]
        return output
    def SumDiagrams(self, amplitudes):
        #This will sum up all the present diagrams in the list of diagrams
        return sum(amplitudes)
    def CountVertices(self, diag):
        #ok I need to write this again
        vertexes=0
        for k in diag.keys():
            if not "out" in str(k):
                vertexes+=1
        return vertexes
    def GetVertices(self, lagrangian):
        #read in the interaction parts of the lagrangian and give a dictionary of vertecies with correponding coupling constants
        #maybe hash down the road to compare models?--This is a far off goal
        vs=dict()
        cc=lagrangian["coupling_constants"]
        for c in cc.keys():
            vs[c]=cc[c]
            #again this is really just reading off a value
        return vs
    def GetPropagators(self, lagrangian):
        #calculate the propagator for each particle that we need to account for
        #return dictionary 
        ps=dict()
        for p in lagrangian["particles"]:
            ps[p]=lagrangian["particles"][p] #the lagrangian data structure should be a dictionary of dictionaries
            #the particle has a mass m, which is the associated propagator
            #this approach will only work for scalars, so further 
        return ps
    def CheckIfVertexValid(self, momenta):
        all_good=False
        weights=[1 for m in momenta]
        s=sum(momenta)
        if s==0:
            all_good=True
        while all_good==False:
            for i in range(len(weights)):
                if momenta[i]==0:
                    continue
                else:
                    weights[i]=-1*weights[i]
                    s=sum([weights[j]*momenta[j] for j in range(len(momenta))])
                    if s==0:
                        all_good=True
                        break
        #this try combinations of inflow or outflow to a vertex
        return all_good

    def AssignValues(self, momenta,verts, frees):
        #this assigns values to the non-free momenta
        max_val=len(momenta)-len(frees)
        vals=list(range(max_val+1))
        m=list(momenta.keys())
        mvs=dict()
        for i in range(len(m)):
            if i in frees:
                fv="free"+str(i)
                momenta[m[i]][1]=fv
                #this has set the free variables
        for v in verts:
            for vl in v:
                if momenta[vl][1]==-1:
                    f_in_v=""
                    needs_random=False
                    for vo in v:
                        if vo==vl:
                            continue
                        else:
                            if "free" in momenta[vo][1]:
                                f_in_v+=momenta[vo][1]
                            elif momenta[vo][1]==-1 :
                                needs_random=True
                            elif momenta[vo][1]==0:
                                continue
                            else:
                                f_in_v+=str(momenta[vo][1])
                    if needs_random:
                        val_to_give=choice(vals)
                        f_in_v+=str(val_to_give)
                    momenta[vl][1]=f_in_v
            return momenta
    def CalculateLoopOrder(self, diagram):
        #takes in a diagram in graph form and outputs the loop order
        l=0
        i=0 
        v=0
        for k in diagram.keys():
            if not "out" in str(k):
                v+=1
                for c in diagram[k][1].keys():
                    if not "out" in str(c):
                        i+=1
        i=i/2
        l=i-v+1
        return l
    def ImposeZeroExternalMomentum(self, diagram):
        #this is the entry point of the CSP to get free parameters to integrate over and external contrivutions
        #this enforces the delta functions that show up at every vertex,
        #nned to change to use jus
        momentum=dict()
        vertices_momenta=list()
        ps=self.propagators
        cc=self.vertices
        all_vertices_are_correct=False
        loop_order=self.CalculateLoopOrder(diagram)
        external_contribution=1
        verts=list()
        for k in diagram.keys():
            if "out" in str(k):
                particle=list(diagram[k][1].values())[0]
                connector=list(diagram[k][1].keys())[0]
                kc=str(k)+":"+str(connector)
                momentum[kc]=[particle, 0]
            #all exterior particles have momentum 0
            #loop order gives the number of free parameters
            #free parameters are given by momentum[n]=[particle, "free"]
            #when checking a vertex with a free parameter, assign "free"=0
            else:
                vert=list()
                external_contribution=external_contribution*-1*cc[diagram[k][0]]
                for nc in diagram[k][1].keys():
                    nkc=str(k)+":"+str(nc)
                    vert.append(nkc)
                    if str(nkc).split(":") in str(momentum.keys()).split(":"):
                        continue
                    else:
                        particle=diagram[k][1][nc]
                        momentum[nkc]=[particle, "free"]
            #Now we have fully setup the dictionary of momenta associated with all lines
            #now we should test assigning values to all but looporder many free partices
                verts.append(vert)
            frees=list()
            if loop_order >0:
                for i in range(int(loop_order)):
                    k=list(momentum.keys())
                    if momentum[k[i]][1]==0:
                        i+=1
                    frees.append(i)

            dmomentum=self.AssignValues(momentum,verts,frees)
            itt=0
            var=-1
        while all_vertices_are_correct==False:
            for k in diagram.keys():
                for nc in diagram[k][1].keys():
                    vertex=list()
                    for m in dmomentum.keys():
                        if str(nc) in str(m) and str(k) in str(m):
                            vm=dmomentum[m][1]
                            if type(vm) is str:
                                if "free" in vm and not "+" in vm:
                                    vm=0
                                elif "free+" in vm:
                                    vm.replace("free+", '')
                                    vm=int(vm)
                                vertex.append(vm)
                            else:
                                continue
                vertices_momenta.append(vertex)
            corrects=list()
            for v in vertices_momenta:
                corrects.append(self.CheckIfVertexValid(v))
            if not all(corrects):
                frees[var]+=itt
                itt+=1
                while itt in frees:
                    itt+=1
                if itt==len(dmomentum):
                    itt=0
                    var+=-1
                k=list(dmomentum.keys())
                while dmomentum[k[itt]][1]==0:
                    itt+=1
                dmomentum=self.AssignValues(momentum, verts, frees)
            else:
                momentum=dmomentum
                all_vertices_are_correct=True
                break
        integrand=list()
        n_var=0
        for m in momentum.keys():
            if "out" in str(m):
                continue
            p=momentum[m][1]
            if p==0:
                continue
            if "free" in str(p):
                if "+" in str(p):
                    p_parts=p.split("+")
                    free_params=list()
                    for k in p_parts:
                        if "free" in str(k):
                            free_params.append(k)
                    add_param=int(p_parts[-1])
                    integrand.append([free_params, add_param, ps[momentum[m][0]]])
                else:
                    n_var+=1
                    var_name="free_"+str(n_var)
                    fp=p.replace("free", var_name)
                    integrand.append([[fp],0, ps[momentum[m][0]]])
            else:
                part=momentum[m][0]
                external_contribution=external_contribution*(-1)/pow(ps[part], 2)
        return [integrand, external_contribution]
    #need to reconstruct for sympy 
    def Integrand(self, xs, xadd, m):
        x=sum(xs)
        return (1/(pow(x+smb.S(xadd),2) +int(pow(smb.S(m),2))))
    def CalculateScatteringAmplitude(self, diagram):
        #this method takes in a diagram in the graph form 
        #then it passes that diagram off to a CSP approach to fix the non-free parameters
        sa=1
        variable_relations, external_contribution=self.ImposeZeroExternalMomentum(diagram)
     #integrand is a list corresponding to an element to integrate for each vertex
     #[free variables depends on, addition to those as just a number, propagator (unused--+)]
        if external_contribution !=0:
            sa=external_contribution
        loop=self.CalculateLoopOrder(diagram)
        #Digest the integrands as they are non-trivial
        #have loop_order many variables to integrate over
        bounds=[[0.000001, self.cutoff]*max(int(loop), 1)]
        #this has given n many copies of the bound
        if loop>0:
            variables=[smb.symbols(x[0][0]) for x in variable_relations if "+" not in str(x[0][0])]
            #Have created an appropriate set of sympy symbols objects
            #loop below assigns each free variable to one of the symbols 
            #then it goes through all of the integands and replaces free variables 
            #this allows for full running of integrand
            for v in variable_relations:
                for st in range(len(v[0])):
                    for vs in variables:
                        if v[0][st]==vs.name:
                            v[0][st]=vs
            integrands=[self.Integrand(v[0], v[1], v[2]) for v in variable_relations]
            pos=-1
            integrand=1
            x=smb.S('x')
            y=smb.S('y')
            expr_multiply=x*y
            for i in range(len(integrands)):
        #        integrand+="*"
                pos=-1 * pos
                integrand=expr_multiply.subs([(x,integrand), (y,smb.sympify(integrands[i]))])
                #print(integrand)
            bounds=[(v, smb.S(0), smb.S(float(self.cutoff))) for v in variables]
            #right now the bounds are coming in as a list, need to be flattened
            #to just read in as (x, 0, cutoff), (y,0, cutoff)...
            #try:
            integral=smb.integrate(pos*smb.sympify(integrand), *bounds)
            #except:
             #   expr=smb.Integral(integrand, *bounds)
              #  print(expr)
               # sys.exit("Bad integral")
            sa=sa*integral.evalf()
           #integrate using sympy then give the proper value of the integral
        #print(sa)
        return sa
    def DrawDiagram(self, diagram):
        #this just takes in the diagram and draws the graph from the dictionary
        #plt.clear()
        Diagram_Graph=nx.Graph()
        particle_labels=dict()
        for k in diagram.keys():
            Diagram_Graph.add_node(k)
        for k in diagram.keys():
            if type(diagram[k][1]) is dict:
                for m in diagram[k][1].keys():
                    Diagram_Graph.add_edge(k,m)
                    particle_labels[(k, m)]=diagram[k][1][m]
        pos=nx.spring_layout(Diagram_Graph)
        #plt.figure()
        nx.draw(Diagram_Graph, pos)
        drawn_graph=nx.draw_networkx_edge_labels(Diagram_Graph, pos,edge_labels=particle_labels)
        plt.axis('off')
        #plt.show(block=False)
        return drawn_graph
    def heuristic(self, diagram, propagator):
        #the heuristic here is given by h=1/SA*(average lambda)/m_prop^2+1-(in+out)/total lines 
        #this last part is a normalized proxy for the loop order
        #all this is saying is really just contributions to scattering amplitude and order of divergence
        #this requires only minimal algebaric calculations 
        #no integrals need to be performed, helping to reduce computational overhead
        avg_constant=0
        total_vertexs=0
        h=1000
        vertexs=self.vertices.copy()
        in_out=0 #this is base level
        if type(diagram) is list:
            for x in diagram[1].keys():
                if "out" in str(x):
                    in_out+=1
            n_lines=1/2*sum([len(diagram[1][x]) for x in diagram[1].keys()])
            in_out_new_contribution=0
            for vk in vertexs.keys():
                if propagator in vk:
                    avg_constant+=vertexs[vk]
                    total_vertexs+=1
                    in_out_new_contribution+=len(vk.split(", "))-1
            if total_vertexs !=0:
                avg_constant=float(avg_constant)/total_vertexs
                in_out_new_contribution=in_out_new_contribution/total_vertexs
            SA=diagram[2] #current scattering amplitude
            m=self.propagators[propagator] #mass of the propagator
            starting_coupling_constant=1
            for node in diagram[1].keys():
                if diagram[1][node][0] != 0:
                    starting_coupling_constant=starting_coupling_constant*vertexs[diagram[1][node][0]]
        #so I need count vertices to return a list of the form 
        #(count of vertexs, product of coupling constants, list of vertexs)
            h=1/SA*(float(starting_coupling_constant)*avg_constant)/pow(m,2)
            n_lines+=in_out_new_contribution
            in_out+=in_out_new_contribution
            h+= in_out/n_lines-1
            return h
        elif type(diagram) is dict:
            for x in diagram.keys():
                if "out" in str(x):
                    in_out+=1
            n_lines=1/2*sum([len(diagram[x]) for x in diagram.keys()])
            in_out_new_contribution=0
            for vk in vertexs.keys():
                if propagator in vk:
                    avg_constant+=vertexs[vk]
                    total_vertexs+=1
                    in_out_new_contribution+=len(vk.split(", "))-1
                if total_vertexs !=0:
                    avg_constant=float(avg_constant)/total_vertexs
                    in_out_new_contribution=in_out_new_contribution/total_vertexs
                SA=self.CalculateScatteringAmplitude(diagram)
                m=self.propagators[propagator] #mass of the propagator
                starting_coupling_constant=1
                for node in diagram.keys():
                    if diagram[node][0] != 0:
                        starting_coupling_constant=starting_coupling_constant*vertexs[diagram[node][0]]
        #so I need count vertices to return a list of the form 
        #(count of vertexs, product of coupling constants, list of vertexs)        
                h=1/SA*(float(starting_coupling_constant)*avg_constant)/pow(m,2)
                n_lines+=in_out_new_contribution
                in_out+=in_out_new_contribution
                h+= in_out/n_lines-1
                return h
            else:
                return 100
        return h
             
        
        
