#This file will holdmath representation python the class to generate feynamn diagrams for searching 
#The generation will need to be a tuple to be hashable
#The hash then should be put in a map with a cost function to improve A* search
from math import comb
from copy import deepcopy
from scipy import integrate
from random import choice
class FeynmanGenerator:
    def __init__(self, lagrangian, cutoff, looporder=1):
        self.order=looporder
        self.cutoff=cutoff
        self.vertices=self.GetVertices(lagrangian)
        self.propagator=self.GetPropagators(lagrangian)
        self.diagrams=self.GenerateDiagrams()
        self.vertex_count=self.CountVertecies()
        self.diagram_output=self.GenerateOutput()
        self.diagram_list=list() #list of the hashs of diagrams
    def UpdateCutoffandVertex(self, new_cutoff, newccs):
        self.cutoff=new_cutoff
        self.vertices=newccs
    def GenerateDiagrams():
        #this will generate diagrams up to loop order given 
        #if loop order=-1, generate until contribution to scattering amplitude is diminneshed by a power of lambda 
        #datastructure for a diagram needs to be 
        #{(string form of incoming to outgoing particles, incoming particles, outgoing particles), set of vertexes, (scattering amplitude numerical, scattering amplitude representational), ingoing-outgoing particles, loop order}
        #this will generate multiple diagrams filling the same set of charracteristics, thus, hash and every incidcdnce of a hash double used, add to the symmetry factor
        scattering_amplitude=0
        diagrams=list()
        hash_diagrams=self.diagram_list()
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
                    outcons=dict()
                    for j in range(len(parts)):
                        out_node="out_"+str(j)
                        out[out_node]=parts[j]
                        out_cons[out_node]=[0, {'1':parts[j]}]
                    out_cons['1']=[v, out]
                    d=tuple([sform, l, rparts], out_cons, vs[v]/sf, 2k-n, 0)
                    #I need to go through all of this and sort out what I need the diagram to have and where
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
    def GenerateNextOrder(self, diagram):
        #need to figure out if I actually want to use this, if so, need to update to match with current Expand
        #what this is doing is actulalay generating the queue
        legs=list()
        children=list()
        for n in diagram[1].keys():
            for m in diagram[1][n].keys():
                leg=str(n)+":"+str(m)
                lgr=str(m)+":"+str(m)
                found=(leg in legs) or (lgr in legs)
                if not found:
                    legs.append(leg)
                    prop_and_nodes=[diagram[1][n][m], n, m]
                    children.append(self.GenerateOutput(diagram, prop_and_nodes))
        return children

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
        ps=self.propagators
        #just take in the graph form of the diagram
        out_diagrams=list()
        external_particles=list()
        diagrams_hash=self.diagram_list
        n_out=0
        for node in diagram.keys():
            if "out" in node:
                n_out+=1
                external_particles.append(list(diagram[node][1].values())[0])
        if second==True:
            external_particles=old_ext
        for v in vs:
            particles_in_vertex=v.split(",")
            if leg in particles_in_vertex:
                new_external=list()
                diagram_to_expand=deepcopy(diagram)
                conections=diagram_to_expand[from_node][1]
                del connections[to_node]
                new_node=len(diagram_to_expand.keys())+1
                connections[new_node]=leg
                diagram_to_expand[from_node][1]=connections
                diagram_to_expand[new_node]=[v,{to_node:leg, from_node:leg}]
                del diagram_to_expand[to_node][1][from_node]
                diagram_to_expand[to_node][1][new_node]=leg
                original_particle=0
                out_n=0
                for r in particles_in_vertex:
                    if original_particle<2 and r==leg:
                        original_particle+=1
                        continue
                    else:
                        out_n+=1
                        node_name="out_"+str(out_n+n_out)
                        diagram_to_expand[new_node][1][node_name]=r
                        diagram_to_expand[node_name]=[0,{new_node: r}]
                        new_external.append(r)
                d=tuple(diagrams_to_expand)
                if not self.Same(new_external, external_particles):
                    if second==False:
                        for m in diagrams_to_expand.keys():
                            for n in diagrams_to_expand[m][1].keys():
                                ds=self.ExpandDiagram(diagrams_to_expand, diagrams_to_expand[m][1][n], m, n, True, external_particles)
                                for dss in ds:
                                    nep=list()
                                    for node in dss.keys():
                                        if "out" in node:
                                            nep.append(list(dss[node][1].values())[0])
                                    if self.Same(nep,external_particles):
                                        tdss=tuple(dss)
                                        if hash(tdss) in diagrams_hash:
                                            continue
                                        else:
                                            out_diagrams.append(dss)
                                            diagrams_hash.append(tdss)
                                    else:
                                        continue

                if hash(d) in diagrams_hash:
                    continue
                else:

                    out_diagrams.append(diagrams_to_expand)
                    diagrams_hash.append(d)
                #so now I've added all the diagrams 
                #now need to put in recombinations
            
                to_recombine=True
                dl=diagrams_to_expand
                while to_recombine:
                    dl=[ self.RecombineExteriors(dx) for dx in dl]
                    dl=[x for y in dl for x in y]
                    if len(dl)==0:
                        to_recombine=False
                        break
                    else:
                        for dd in dl:
                            d=tuple(dd)
                            new_external=list()
                            for node in dd:
                                if "out" in node:
                                    new_external.append(list(dd[node][1].values())[0])
                            if not self.Same(new_external, external_particles):
                                continue
                            if hash(d) in diagrams_hash:
                                continue
                            else:
                                out_diagrams.append(dd)
                                diagrams_hash.append(d)
                return out_diagrams


    def RecombineExteriors(self, diagram):
        #recombines exterior propagators
        k=[x for x in diagram.keys() if "out" in x]
        found_recomb=False
        diagrams_out=list()
        dummy_diagram=deepcopy(diagram)
        d_hash=list()
        for i in range(len(k)):
            for j in range(len(k)):
                if i==j: 
                    continue
                else:
                    if  list(dummy_diagram[k[i]].values())[0]==list(dummy_diagram[k[j]].values())[0]:
                        found_recomb=True
                        c=list(dummy_diagram[k[j]].keys())[0]
                        b=list(dummy_diagram[k[i]].keys())[0]
                        p=list(dummy_diagram[k[i]].values())[0]
                        dummy_diagram[b][c]=p
                        dummy_diagram[c][b]=p
                        del dummy_diagram[k[j]]
                        del dummy_diagram[k[i]]
                        
                        d=tuple(dummy_diagram)
                        if hash(d) in d_hash:
                            continue
                        else:
                            diagrams_out.append(dummy_diagram)
                            d_hash.append(d)
        return diagrams_out


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
            if not "out" in k:
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
    def GetPropagators(self, lagrangian)
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
        while all_good=False
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
        vals=[0:max_val]
        m=list(momenta.keys())
        mvs=dict()
        for i in range(len(m)):
            if i in frees:
                fv="free"+str(i)
                momenta[m[i]][1]=fv
                #this has set the free variables
        for v in verts:
            for vl in v
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
                            else
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
            if not "out" in k:
                v+=1
                for c in diagram[k][1].keys():
                    if not "out" in c:
                        i+=1
        i=i/2
        l=i-v+1
        return l
    def ImposeZeroExternalMomentum(self, diagram):
        #this is the entry point of the CSP to get free parameters to integrate over and external contrivutions
        #this enforces the delta functions that show up at every vertex
        momentum=dict()
        vertices_momenta=list()
        ps=self.propagators
        cc=self.vertices
        all_vertices_are_correct=False
        loop_order=self.CalculateLoopOrder(diagram)
        external_cotnribution=1
        verts=list()
        for k in diagram.keys():
            if "out" in k:
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
                external_contribution=external_contribution*-1*cc[k][0]
                for nc in diagram[k]:
                    nkc=str(k)+":"+str(nc)
                    vert.append(nkc)
                    if nkc.split(":") in momentum.keys().split(":"):
                        continue
                    else:
                        particle=diagram[k][1][nc]
                        momentum[nkc]=[particle, "free"]
            #Now we have fully setup the dictionary of momenta associated with all lines
            #now we should test assigning values to all but looporder many free partices
                verts.append(vert)
            frees=list()
            f=len(momentum) % loop_order
            for i in range(len(loop_order)):
                k=list(momentum.keys())
                if momentum[k[i]][1]==0:
                    i+=1
                frees.append(i)

            dmomentum=self.AssignValues(momentum,verts,frees)
            itt=0
            var=-1
        while all_vertices_are_correct==False:
            for k in diagrams.keys():
                for nc in diagrams[k].keys():
                    vertex=list()
                    for m in dmomentum.keys()
                        if nc in m and k in m:
                            vm=dmomentum[m][1]
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
        for m in momentum.keys():
            if "out" in m:
                continue
            p=momentum[m][1]
            if "free" in p:
                if "+" in p:
                    p_parts=p.split("+")
                    free_params=list()
                    for k in p_parts:
                        if "free" in k:
                            free_params.append(k)
                    add_param=int(p_parts[-1])
                    integrand.append([free_params, add_param, ps[momementum[m][0]]])
                else:
                    fp=p.replace("free", '')
                    integrand.append([[fp],0, ps[momentum[m][0]])
            else:
                part=momentum[m][0]
                external_contribution=external_contribution*(-1)/pow(ps[part], 2)
        return [integrand, external_contribution]
    def Integrand(self, xs, x_add, m):
        x=sum(xs)
        return -1*/(pow(x +x_add,2) +pow(m,2))
    def IntegrandProduct(self, integrands):
        I=1
        for i in integrands:
            I=I*self.Integrand(i[0], i[1], i[2])
        return I
    def CalculateScatteringAmplitude(self, diagram):
        #this method takes in a diagram in the graph form 
        #then it passes that diagram off to a CSP approach to fix the non-free parameters
        sa=1
        pi=3.14159
        integrands, external_contribution=self.ImposeZeroExternalMomentum(diagram)
        sa=external_contrbution
        loop=self.CalculateLoopOrder(diagram)
        #now I need to digest the integrands as they are non-trivial
        #have loop_order many variables to integrate over
        bounds=[[0, self.cutoff]]*loop
        #this has given n many copies of the bound
        sa=sa*integrate.nquad(self.IntegrandProduct, bounds, args=(*integrands))    
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
        

             
        
        
