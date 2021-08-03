import xlrd
import codecs
from collections import defaultdict
def get_disease_mesh_cni(target_disease):
    data = xlrd.open_workbook(target_disease)
    data.sheet_names()
    table = data.sheet_by_name('target_cancer')
    rows = table.nrows
    cols = table.ncols
    hasdir = defaultdict(dir)
    for rownum in range(1,rows):
        row_values = [item.strip() for item in table.row_values(rownum)]
        diseasename = row_values[0]
        mesh = row_values[1]
        cni = [item for item in row_values[2:] if item!="-" and item!='']
        hasdir[diseasename]={mesh:cni}
    return hasdir

def mesh_miro(micro_file,disease_mesh,output):
    wf = codecs.open(output,'w',encoding='utf-8')
    wf.write("DiseseID\tDirection\tMicroID\n")
    with open(micro_file,'r') as rf:
        rf.readline()
        for line in rf:
            contents = line.split('\t')
            micro = contents[0]
            disease = contents[1]
            direction = contents[2]
            mesh_dise = [item for item in disease_mesh[disease].keys() if item!='-']
            if direction=="up":
                direction=1
            elif direction=="down":
                direction=-1
            else:
                print(direction)
                direction=None
            if len(mesh_dise)!=0 and direction!=None:
                wf.write("{}\t{}\t{}\n".format("MESH:"+mesh_dise[0],direction,micro))
    wf.close()
 
def main():
    target_disease = "./target_disease.xlsx"
    micro_file = "./microrna/miRCancerOctober2019.txt"
    output = "./microrna/mesh2mcro.txt"
    
    disease_mesh = get_disease_mesh_cni(target_disease)
    disease_mesh = {key:value for key,value in disease_mesh.items()}
    mesh_miro(micro_file,disease_mesh,output)

if __name__=="__main__":
    main()