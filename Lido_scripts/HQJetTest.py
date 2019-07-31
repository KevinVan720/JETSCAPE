import subprocess
from subprocess import Popen
pTbins=[5, 15, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 220, 250]
#pTbins=[100,110]
#for i in range(len(pTbins)-1):
#    subprocess.run(["./HQJetTest", str(pTbins[i]), str(pTbins[i+1])])

#commands = ["./HQJetTest " + str(pTbins[i])+" "+str(pTbins[i+1]) for i in range(len(pTbins)-1)]
procs = [ Popen(["./HQJetTest", str(pTbins[i]), str(pTbins[i+1]), str(200000)]) for i in range(len(pTbins)-1) ]
for p in procs:
   p.wait()