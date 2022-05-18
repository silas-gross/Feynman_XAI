#This file is to perform the main searching over the Feynman diagrams, LAO* approach is goal
#CSP is also to be put-into this search, fill back constraints to cut methods
from Feynman_generator  import FeynmanGenerator
from Renormalization    import Willson
from Renormalization    import Sythetic
class FeynmanSearch:
    def __init__(self, Lagrangian):
        cutoff_initial=0
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
            while good_change=False:
                dcc=-self.cutoff*deltasa/dco
                good_change=self.SynthConstraint(self.cutoff, dco, dcc, vert, deltasa)
                if good_change:
                    break
                else:
                    dco+=-0.05*dco
            self.ccs[verts]+=dcc
            self.cutoff+=-dco


        #Willsonian here would be beta function which imposes that 1/cc-1/new cc=change in SA*log(1/cutoff)
        if method="WR":

        
    def PerformSearchWillson(self, diagram):
#this right now needs to be fleshed out more, and I need to figure out what a sucessful goal can be defined as 
    def PerformSearchRunning(self, diagram)
    def DecidePointFunction(self):
        generator=self.diagram_base
        diagrams=generator.diagrams
        correction=dict()
        for d in diagrams:
            cs=generator.GenerateNextOrder(d[1])
            c=sum(cs)/len(cs)
            dh=hash(d)
            correction[dh]=c/d[2]
        diagram_to_use=diagrams[0]
        for d in diagrams:
            if correction[hash(d)]<=correction[hash(diagram_to_use)]:
                diagram_to_use=d
            else:
                continue
        return diagram_to_use
