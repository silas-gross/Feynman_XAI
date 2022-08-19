#This file is to perform the main searching over the Feynman diagrams, LAO* approach is goal
#CSP is also to be put-into this search, fill back constraints to cut methods
from Feynman_generator  import FeynmanGenerator
#from Renormalization    import Willson
#from Renormalization    import Sythetic
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
        self.cutorr=4*ci
    def SynthConstraint(self, co, dco, dcc, vert, deltasa):
        lhs=co*deltasa
        rhs=(co-dco)*sc*(self.ccs[vert]+dco)
        if abs(lhs-rhs) < lhs/abs(dco-dcc):
            return True
        else:
            return False
    def ImposeConstraints(self, method, verts, deltasa, diagram):
        #for SC/RC the constriaint is that 
        #new cc = cc -change in SA *cuttoff/change in cutoff
        if method=="SC":
            dco=0.99*self.cutoff*deltasa/self.ccs[verts]
            dcc=0
            good_change=False
            while good_change==False:
                dcc=-self.cutoff*deltasa/dco
                good_change=self.SynthConstraint(self.cutoff, dco, dcc, vert, deltasa)
                if good_change:
                    break
                else:
                    dco+=-0.05*dco
            if self.ccs[verts]+dcc<=0 or self.cutoff-dco<=0:
                return False
            else:
                self.ccs[verts]+=dcc
                self.cutoff+=-dco
                return True



        #Willsonian here would be beta function which imposes that 1/cc-1/new cc=change in SA*log(1/cutoff)
        if method=="WR":
            co=self.cutoff
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
        generator.DrawDiagram(di)
        plt.show()
        at_goal=False
        children=generator.GenerateNextOrder(d)
        queue=dict()
        scattering_amp=generator.CalculateScatteringAmplitude(d)
        for c in children:
            if c[2] in queue.keys():
                queue[c[2]].append([c[0], c[1]])
            else:
                queue[c[2]]=[[c[0], c[1]]]
        hs=list(queue.keys())
        highest_priority=min(hs)
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
                knew=d1.keys()
                for k in knew:
                    if not k in kold:
                        good=self.ImposeConstraints(method, k, deltasa[-1],d1)
                        if not good:
                            break
                        else:
                            generator.UpdateCutoffandVertex(self.cutoff, self.css)
            sa=sum(deltasa)
            scattering_amp+=sa
            if sa<highest_priority:
                at_goal=True
                return [cd[0], scattering_amplitude]
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
            cs=generator.GenerateNextOrder(d[1])
            c=sum([x[2] for x in cs])/len(cs)
            dh=hash(str(d))
            correction[dh]=c/d[2]
        diagram_to_use=diagrams[0]
        for d in diagrams:
            if correction[hash(str(d))]<=correction[hash(str(diagram_to_use))]:
                diagram_to_use=d
            else:
                continue
        self.diagram_to_use=diagram_to_use
        generator.DrawDiagram(diagram_to_use[1])
        plt.show()
        return diagram_to_use
