import iFeatureOmegaCLI
import os
import glob

folder = glob.glob("combo/*.fasta")
for fl in folder:
        print(fl)
        prot = iFeatureOmegaCLI.iProtein(fl)
        prot.import_parameters('protein_params.json')
        prot.get_descriptor('DPC type 1')
        prot.to_csv(fl+'_2mer.csv', index=True, header=True)
        prot.get_descriptor('PAAC')
        prot.to_csv(fl+'_pseaac.csv', index=True, header=True)
        prot.get_descriptor('CTriad')
        prot.to_csv(fl+'_ctriad.csv', index=True, header=True)
        prot.get_descriptor('CTDC')
        prot.to_csv(fl+'_ctdc.csv', index=True, header=True)
        prot.get_descriptor('CTDT')
        prot.to_csv(fl+'_ctdt.csv', index=True, header=True)
        prot.get_descriptor('CTDD')
        prot.to_csv(fl+'_ctdd.csv', index=True, header=True)

print("FINITO!")

## comment added
