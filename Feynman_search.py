#This file is to perform the main searching over the Feynman diagrams, A* approach to search for most productive branches
#CSP is also to be put-into this search, fill back constraints to cut methods
from Feynman_generator  import FeynmanGenerator
import math
import matplotlib.pyplot as plt

class FeynmanSearch:
    def __init__(self, Lagrangian):
        cutoff_initial=0
        self.l=Lagrangian
        self.diagram_to_use=list()
        for m in Lagrangian["particles"].values():
            cutoff_initial+=m*m
        cutoff_intial=4*cutoff_initial
        self.cutoff=cutoff_initial
        self.ccs=Lagrangian["coupling_constants"]
        self.masses=Lagrangian["particles"]
        #note that this is the sum of the square of the energies of a kinetic energy=m
        #so this is a intial very energetic system
        self.diagram_base=FeynmanGenerator(Lagrangian, cutoff_intial)
        #goal is delta cc <heuristic cost of lowest
    def ResetValues(self):
        self.css=self.l["coupling_constants"]
        self.masses=self.l["particles"].values()
        ci=sum([m*m for m in self.masses])
        self.cutoff=4*ci
    def SynthConstraint(self, co, dco, dcc, vert, deltasa):
        #this function makes sure that all of the changes stay manifestly positive
        lhs=co/dco
        rhs=self.ccs[vert]/deltasa
        #print(rhs)
        if lhs-rhs < lhs:
            return True
        else:
            return False
    def ImposeConstraints(self, method, verts, deltasa, diagram):
        #for SC/RC the constriaint is that 
        #new cc = cc -change in SA *cuttoff/change in cutoff
        #needs a mass change as well I think
        if method=="SC":
            scale=0.01 
            #this is introduced to allow for better scaning over phase space
            #prevents runaway checks of nonsensicial changes to coulping constants
            if verts not in self.ccs.keys():
                return False
            dco=0.99*self.cutoff
            print("initial change to cutoff: ", dco)
            dcc=0
            good_change=False
            while good_change==False:
                dcc=-1*scale*self.cutoff*deltasa/(dco)
                good_change=self.SynthConstraint(self.cutoff, dco, dcc, verts, deltasa)
                #print(self.cutoff)
                if good_change:
                    break
                else:
                    dco+=-0.01*dco
                    if abs(dcc) >abs(self.ccs[verts]):
                        scale= self.ccs[verts]/dcc
                     #   break
                    #    print("bad change", dco,dcc)
            if self.ccs[verts]+dcc<=0 or self.cutoff-dco<=0:
                print("bad change, issue is verts? ", self.ccs[verts]+dcc, "or cutoff?", self.cutoff-dco)
                return False
            else:
                self.ccs[verts]+=dcc
                self.cutoff+=-dco
                print("changed cuttoff: ", self.cutoff)
                return True



        #Willsonian here would be beta function which imposes that 1/cc-1/new cc=change in SA*log(1/cutoff)
        if method=="WR":
            co=self.cutoff
            print(co)
            vs=list()
            for d in diagram.keys():
                vs.append(diagram[d][0])
            for v in vs:
                if not v in self.css.keys():
                    continue
                c=1/(1/self.ccs[v] +deltasa*math.log(co))
                if c<=0 or c>=self.ccs[v]:
                    return False
                else:
                    self.ccs[v]=c
            return True

    def PerformSearch(self, method):
        #get priority queue from Feynman Generator
        #first reset the coupling constants to the original values
        self.ResetValues()
        generator=self.diagram_base
        di=self.diagram_to_use
        d=di[1]
        generator.DrawDiagram(d)
        #plt.show()
        children=generator.GenerateNextOrder(d)
        queue=dict()
        scattering_amp=generator.CalculateScatteringAmplitude(d)
        for c in children:
            if c[2] in queue.keys():
                queue[c[2]].append([c[0], c[1]])
            else:
                queue[c[2]]=[[c[0], c[1]]]
        hs=list(queue.keys())
        htemp=None
        highest_priority=0
        for h in hs:
            ha=abs(h)
            if htemp is None or ha<htemp:
                highest_priority=h
                htemp=ha
        #print(highest_priority)
        while len(queue)>0:
            if len(queue[highest_priority])==0:
                queue.pop(highest_priority)
                hs=list(queue.keys())
                if len(hs)==0:
                    break
                hs.sort(reverse=True)
                highest_priority=hs[0]
            cd=queue[highest_priority].pop()
            diags=generator.ExpandDiagram(cd[0], *cd[1])
            kold=cd[0].keys()
            deltasa=[]
            for d1 in diags:
                deltasa.append(generator.CalculateScatteringAmplitude(d1))
                #print("diagram change is ", deltasa[-1])
                knew=d1.keys()
                #print(knew)
                for k in knew:
                    if not k in kold:
                        good=self.ImposeConstraints(method, d1[k][0], deltasa[-1],d1)
                        if not good:
                            #print("going to do a break", len(queue))
                            break
                        else:
                            print(self.cutoff)
                            generator.UpdateCutoffandVertex(self.cutoff, self.css)
                    else:
                        continue
            sa=sum(deltasa)
            scattering_amp+=sa
            #print(highest_priority)
            if sa<abs(highest_priority):
                return [cd[0], scattering_amp]
            else:
                child=generator.GenerateNextOrder(cd[0])
                for c in child:
                    if c[2] in queue.keys():
                        queue[c[2]].append([c[0], c[1]])
                    else:
                        queue[c[2]]=[[c[0], c[1]]]
                        highest_priority=max(c[2], highest_priority)

                

    def DecidePointFunction(self):
        generator=self.diagram_base
        diagrams=generator.diagrams
        correction=dict()
        for d in diagrams:
            #print(d)
            cs=generator.GenerateNextOrder(d)
            c=sum([x[2] for x in cs])/len(cs)
            dh=hash(str(d))
            correction[dh]=c/d[2]
        diagram_to_use=diagrams[0]
        for d in diagrams:
            if correction[hash(str(d))]>=correction[hash(str(diagram_to_use))]:
                diagram_to_use=d
            else:
                continue
        self.diagram_to_use=diagram_to_use
        generator.DrawDiagram(diagram_to_use[1])
        #plt.show()
        return diagram_to_use
