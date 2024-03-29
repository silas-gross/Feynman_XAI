#this is the enerty point and contains the GUI
import tk as tk
from Feynman_search import FeynmanSearch
import time 
import copy
import multiprocessing as mp

#top=tk.Tk()
def SubmitLagrangian():
    pass
def EnterLagrangian():
    w=tk.Canvas(top, cursor="dot", height=600, width=500)
    w.create_text(50, 250, text="Enter Lagrangian", font=('calibre', 15, 'bold'))
    w.create_text(50, 300, text="Number of Real Scalar Fields: ", font=('calibre', 15))
    w.create_text(50, 350, text="Number of Complex Scalar Fields: ", font=('calibre', 15))
    w.pack()
    w.Draw()

#
#def DrawFeynmanDiagram(diag, title=''):

def LagrangianInput():
    #command line input for Lagrangian
    #idealy want to move this over to GUI but will see if I get the time
    nr=0
    nc=0
    mr=[]
    mc=[]
    ni=0
    ccs=dict()
    #print("Enter the Lagrangian (real and complex scalars only): \n")
    #print("Give number of real scalar fields: ")
    #nr=int(input())
    nr=1
    #print("\n Give number of complex scalar fields: ")
    #nc=int(input())
    nc=1
    mr=[12.5]
    #for i in range(nr):
     #   print("\n Give mass of real scalar field %d in GeV: ", (i+1))
      #  mr.append(float(input()))
    #for i in range(nc):
     #   print("\n Give mass of complex scalar field %d in GeV: ", (i+1))
     #   mc.append(float(input()))
    mc=[125.6]
    #print("\n Give number of interaction terms in Lagrangian: ")
    #ni=int(input())
    ni=2
    #print("\n Give interaction terms with particles comma seperated. Use phi_i for real scalars, psi_i for complex")
    #print("\n For example a vertex with two real scalar fields and two particles of the same complex scalar would be \" phi_1, phi_2, psi_1, psi_2 \" ")
    #for i in range(ni):
    #    print("\n Give particles in interaction term %d seperated by commas.", (i+1))
    #    v=str(input())
    #    print("\n Give value of coupling constant: ")
    #    cc=float(input())
     #   ccs[v]=cc
    v_1="phi_1, phi_1, phi_1, phi_1"
    cc_1=125
    v_2="phi_1, psi_1, psi_1"
    cc_2=1000
    v_3="phi_1, phi_1, psi_1"
    cc_3=12
    #ccs[v_1]=cc_1
    ccs[v_2]=cc_2
    #ccs[v_3]=cc_3
    particles=dict()
    for i in range(nr):
        pname="phi_"+str(i+1)
        particles[pname]=mr[i]
    for i in range(nc):
        pname="psi_"+str(i+1)
        particles[pname]=mc[i]
    lagrangian={"particles":particles, "coupling_constants":ccs}
    return lagrangian
#w=tk.Canvas(top, cursor="dot", height=400, width=500)
#w.create_text(100, 250, text="Hello, welcome to renormalization!")
#w.create_text(150, 250, text="This program takes in a Lagrangian, searches for renormalization by expanding along lines in n-point functions of feynman diagrams \n Right now this is limited to scalar field theories, but will be modified to later allow for spinnors")
#w.create_text(250, 250, text="Then renormalization will be carried through using Willson's approach and synthetic coupling, showing equivielnce \n This equivilence is useful in connecting back to xAI and learning (see report)")
#w.pack()
#n=tk.Button(top, text="Move onto Enter Lagrangian", command=EnterLagrangian(w))
#n.pack()
#top.mainloop()

#Can't get this working properly right now so that is not a priority

print("Hello, welcome to renormalization! \n This program takes in a Lagrangian, searches for renormalization by expanding along lines in n-point functions of feynman diagrams \n Right now this is limited to scalar field theories, but will be modified to later allow for spinnors \n Then renormalization will be carried through using Willson's approach and synthetic coupling, showing equivielnce \n This equivilence is useful in connecting back to xAI and learning (see report) \n \n \n")
print("\n To get started please follow the instructions to enter the Lagrangian")
L=LagrangianInput()
L_fixed=copy.deepcopy(L)
print("\n Thank you for entering lagrangian, now we will expand all the vertex and tree level diagrams")
print("\n This expansion chooses the n-point function that we will investigate for corrections")
sc=FeynmanSearch(L)
ic={c:L["coupling_constants"][c] for c in L["coupling_constants"]}
vertexdiagrams=sc.diagram_base.diagrams
graphdiagrams=dict()
visual_diagrams=list()
for vd in vertexdiagrams:
    g=vd[1]
    v=""
    for m in g.keys():
        if not "out" in m:
            v=g[m][0]
    graphdiagrams[v]=g
for d in graphdiagrams.keys():
    visual_diagrams.append([graphdiagrams[d], str(d)])
    #this gives the diagram in a visulaizable form, can then call the Draw method
    #from Feynman generator
    
diagram=sc.DecidePointFunction()
print("\n The expansion will happen on the diagrams corresponding to traveling between states: " + diagram[0][0])
[wilson_answer_diagram, wSA]=sc.PerformSearch("WR")
willson_couplings=sc.css
willson_masses=sc.l["particles"]
wco=sc.cutoff
sc.l=copy.deepcopy(L_fixed)
sc.ResetValues()
[sc_answer_diagram, scSA]=sc.PerformSearch("SC")
sc_couplings=sc.css
sc_masses=sc.l["particles"]
sc_cutoff=sc.cutoff
print("\n finished running, results have been output to file output.txt")
outfile=open("output.txt", "w")
outfile.write("Willson Renormalization answer \n \n")
outfile.write("diagram graph: " +str(wilson_answer_diagram))
outfile.write("\n initial coupling: "+str(ic))
outfile.write("\n final scattering amplitude: " +str(wSA))
outfile.write("\n vertex couplings: " +str(willson_couplings))
outfile.write("\n masses: " +str(willson_masses))
outfile.write("\n cutoff momenta: "+ str(wco))
outfile.write("\n \n Synthetic Coupling Renormalization answer \n \n")
outfile.write("diagram graph: " +str(sc_answer_diagram))
outfile.write("\n final scattering amplitude: " +str(scSA))
outfile.write("\n vertex couplings: " +str(sc_couplings))
outfile.write("\n cutoff momenta: " +str(sc_cutoff))
outfile.write("\n masses: "+str(sc_masses))
outfile.close()
print("\n Bye!")

