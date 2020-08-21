import sys

if len(sys.argv) != 5:
    print("Usage: python smk2wdl.py ProjectID SubProjectID smk.sample.txt wdl.sample.txt")
else:
    projectID, subprojectID, smk, wdl = sys.argv[1:]

f_o = open(wdl, "w")
f_o.write("#ProjectID\tSubProjectID\tSampleID\tLibraryID\tFQ1\tFQ2\tMark\tInsertSize\tSubSampleID\n")

for i in open(smk).readlines()[1:]:
    item = i.strip().split("\t")
    nstr = "%s\t%s\t%s\t%s\t%s\t%s\t1\t500\tmeta\n" %(projectID, subprojectID, item[0], item[0], item[1], item[2])
    f_o.write(nstr)

f_o.close()
