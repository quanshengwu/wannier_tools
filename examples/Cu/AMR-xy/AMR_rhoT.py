# Merge data of the same temperature from different folders.
import os
import re

def extract_different_temp_data(f, ca_num, output_dir):
    """
    :param f:
    :param ca_num
    :return:
    """
    with open(f, 'r', encoding='utf-8') as data_file:
        line = data_file.readline()
        while line:
            if '#  T' in line:
                line=line.strip()
                tmp = line.split(' ')[-2]
                output_file = './{}/{}K.dat'.format(output_dir, tmp)
                with open(output_file, 'a+', encoding='utf-8') as out:
                    out.write('# {}\n'.format(ca_num))
                    line = data_file.readline()
                    while line:
                        if '#  T' not in line:
                            out.write(line)
                            line = data_file.readline()
                        else:
                            break
            else:
                line = data_file.readline()


if __name__ == '__main__':
    number=1 #Btheta as the variable, number=0; Bphi as the variable, number=1
    dirs = os.listdir('.')
    if not os.path.exists('./rho'):
        os.mkdir('./rho')
    files= os.listdir('./rho')
    for file in files:
        if file.endswith('K.dat'):
            file_path = os.path.join('./rho', file)
            os.remove(file_path)

    tmp_list = []
    for dir_name in dirs:
        if 'Btheta' in dir_name:
            angle=re.findall(r'\d+',dir_name)
            tmp=int(angle[number])
            tmp_list.append(tmp)
    tmp_list.sort()
    for dir_name_num in tmp_list:
        ca_num = str(int(dir_name_num))
        #dir_name = 'Btheta' + str(dir_name_num) + 'Bphi90'
        dir_name = 'Btheta90'+ 'Bphi' + str(dir_name_num) 
        extract_different_temp_data('./{}/rho_total_mu_0.00eV.dat'.format(dir_name), ca_num=ca_num, output_dir='rho')
