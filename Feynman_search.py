#This file is to perform the main searching over the Feynman diagrams, A* approach to search for most productive branches
#CSP is also to be put-into this search, fill back constraints to cut methods
from Feynman_generator  import FeynmanGenerator
import math
import matplotlib.pyplot as plt
import multiprocessing as mt
import copy 
import time 

class PlottingThread(): #a seperate thread to do plotting of the diagrams to not hold up exectution 
    def __init__(self, generator, diagram, name='plt_thread'):
        self.diagram=diagram
        self.gen=generator
#        super(PlottingThread, self).__init__(name=name)
#        self.start()
    def run(self):
        self.gen.DrawDiagram(self.diagram)
        plt.show()
    def UpdateDiagram(self, d):
        plt.clf()
        self.diagram=d
        self.gen.DrawDiagram(self.diagram)
        plt.show()
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
        rhs=(self.ccs[vert]-dcc)/deltasa
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
            #print("initial change to cutoff: ", dco)
            dcc=0
            good_change=False
            while good_change==False:
                dcc=-1*self.cutoff*deltasa/(dco)
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
            if self.cutoff-dco<=0:#self.ccs[verts]+2*dcc<=0 or self.cutoff-dco<=0:
                print("bad change, issue is verts? ", self.ccs[verts]+dcc, "or cutoff?", self.cutoff-dco)
                return False
            else:
                self.ccs[verts]+=dcc
                self.cutoff+=-dco
#                print("changed cuttoff: ", self.cutoff)
                return True



        #Willsonian here would be beta function which imposes that 1/cc-1/new cc=change in SA*log(1/cutoff)
        if method=="WR":
            co=self.cutoff
            #print(co)
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
    def Calculate(self, diagram, generator, scattering_amplitude):
        scattering_amplitude.append(generator.CalculateScatteringAmplitude(diagram))

    def PerformSearch(self, method):
        #get priority queue from Feynman Generator
        #first reset the coupling constants to the original values
        self.ResetValues()
        t_start=time.time()
        generator=self.diagram_base
        di=self.diagram_to_use
        d=di[1]
        old_diagr=copy.deepcopy(di)
        kthread=mt.Process(target=PlottingThread, args=(generator, d))
        #generator.DrawDiagram(d)
        #plt.show()
        kthread.start()
        children=generator.GenerateNextOrder(d)
#        print(children)
        queue=dict()
        scattering_amp=generator.CalculateScatteringAmplitude(d)
        for c in children:
            if c[2] in queue.keys():
                queue[c[2]].append([c[0], c[1]])
            else:
                queue[c[2]]=[[c[0], c[1]]]
        #print(queue.keys())
        hs=list(queue.keys())
        htemp=None
        highest_priority=0
        for h in hs:
            ha=abs(h)
            if htemp is None or ha<htemp:
                highest_priority=h
                htemp=ha
        print("length of queue: "+str(len(queue)))
        while len(queue)>0:
            if time.time() >= t_start+900:
                #timeout if we are beyond 30 minues
                while len(queue[highest_priority])==0:
                    queue.pop(highest_priority)
                break
            if len(queue[highest_priority])==0:
                queue.pop(highest_priority)
                hs=list(queue.keys())
                if len(hs)==0:
                    break
                hs.sort(reverse=True)
                highest_priority=hs[0]
            cd=queue[highest_priority].pop()
            print("number of verticies in top diagram is " + str(generator.CountVertices(cd[0])))
            diags=generator.ExpandDiagram(cd[0], *cd[1])
        #    kthread.UpdateDiagram(cd[0])
            kold=cd[0].keys()
            deltasa=[]
            for d1 in diags:
                if d1==old_diagr:
                    break
                if generator.heuristic(d1, "psi_1") < 0.2 * scattering_amp:
                    deltasa.append(0.01)
                else:
                    p=mt.Process(target=self.Calculate, args=(d1, generator, deltasa))
                    p.start()
                    p.join(30) #wait one minuuee to attempt calculation of scattering amplitude
                    if p.is_alive():
                        p.kill()
                        p.join()
                        deltasa.append(0.01)
                    else:
                        deltasa.append(generator.CalculateScatteringAmplitude(d1))
                knew=d1.keys()
                #print(knew)
                for k in knew:
                    if not k in kold:
                        good=self.ImposeConstraints(method, d1[k][0], deltasa[-1],d1)
                        if not good:
                            #print("going to do a break", len(queue))
                            break
                        else:
                       #     print(self.cutoff)
                            generator.UpdateCutoffandVertex(self.cutoff, self.css)
                    else:
                        continue
#            print("Outside of the loop now")
            sa=sum(deltasa)
            scattering_amp+=sa
 #           print(highest_priority)
            if abs(sa)<0.9*abs(highest_priority):
                kthread.UpdateDiagram(cd[0])
                kthread.join()
                generator.DrawDiagram(cd[0])
                plt.show()
                return [cd[0], scattering_amp]
                break
            else:
                children=generator.ExpandDiagram(cd[0], *cd[1])
                d=cd[0]
                old_diags=copy.deepcopy(d)
                for child_a in children:
                    child=generator.GenerateNextOrder(child_a)
                    for c in child:
                        #generator.DrawDiagram(c)
                        #plt.show()
                        print("Childer have this number of vertices: "+ str(generator.CountVertices(c[0])))
                        if hash(str(c)) in generator.diagram_list: continue
                        if generator.CountVertices(c[0]) == 1: continue
                        if c[2] in queue.keys():
                            queue[c[2]].append([c[0], c[1]])
                        else:
                            queue[c[2]]=[[c[0], c[1]]]
                            highest_priority=max(c[2], highest_priority)
                print("queue length is now " + str(len(queue)))
        generator.DrawDiagram(d)
        plt.show()
        return[d, scattering_amp]        

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
        plt.show()
        generator.exterior=diagram_to_use[0][1]+diagram_to_use[0][2]
        return diagram_to_use
