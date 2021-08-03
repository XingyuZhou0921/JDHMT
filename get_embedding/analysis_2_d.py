# encoding:utf-8

import os
def split_out_disease(output,inputdata,meshdisease,xmin,xmax,ymin,ymax,index):
    basename = os.path.basename(inputdata)+str(index)
    wf = open(os.path.join(output,basename),'w')
    with open(inputdata) as rf:
        rf.readline()
        for line in rf:
            contents = line.strip().split(",")
            if contents[0]=="0":
                x=float(contents[2])
                y = float(contents[3])
                if x>=xmin and x<=xmax and y>=ymin and y<=ymax:
                    mesh = contents[1]
                    wf.write("{}\t{}\n".format(mesh,meshdisease[mesh]))

def mapmesh2diseasename(orgindata):
    mesh2diseasename = {}
    with open(orgindata) as rf:
        for line in rf:
            contents = line.strip().split("\t")
            mesh = "MESH:"+contents[1]
            name = contents[0]
            mesh2diseasename[mesh]=name
    return mesh2diseasename

        
        


def main():
    output = "/mnt/disk1/RESCALSIDE/RESCAL_SIDEINFORMATION/t-SNE_analysis/"
    inputdata = "/mnt/disk1/RESCALSIDE/RESCAL_SIDEINFORMATION/t_SNE/single_rescal_unfilter.csv"
    diseas_mesh = "/mnt/disk1/RESCALSIDE/DataForWA/diseasename2mesh.txt"
    targetdisease = "/mnt/disk1/RESCALSIDE/RESCAL_SIDEINFORMATION/targetdisease.txt"
    xmin = 0
    xmax = 20
    ymin = 20
    ymax = 40
    meshdisease1 = mapmesh2diseasename(diseas_mesh)
    meshdisease2 = mapmesh2diseasename(targetdisease)
    meshdisease = {**meshdisease1,**meshdisease2}
    split_out_disease(output,inputdata,meshdisease,xmin,xmax,ymin,ymax,index=2)

if __name__=="__main__":
    main()


