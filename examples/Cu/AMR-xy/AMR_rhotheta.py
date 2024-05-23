import os
import re
import numpy as np
def extract_rho_theta(f,dtheta,num_show,interval_line,choose,output_file, output_dir):
    index = -1
    for _ in range(num_show):
       theta=-dtheta
       index = index+1
       with open(f, 'r', encoding='utf-8') as data_file:
          line = data_file.readline()
          outfile = './{}/{}'.format(output_dir,output_file)
          with open(outfile, 'a+', encoding='utf-8') as out:
               out.write('#Btau = {:.2f}\n'.format(choose[index]))
               out.write('\n')
               while line:
                  if 'BTau' in line:
                     theta=theta+dtheta
                     for __ in range(index*interval_line+1):
                         line=data_file.readline()
                     parts=line.split()
                     formparts=['{:>14}'.format(part) for part in parts]
                     formparts[0]='{:>5}'.format(str(theta))
                     line='   '.join(formparts)+'\n'
                     out.write(line)
                     line = data_file.readline()
                  else:
                     line = data_file.readline()


if __name__ == '__main__':
    theta_interval=15
    Btau_show=6
    Btau_num=101
    Btau_max=10
    Btau_list=np.linspace(0,Btau_max,Btau_num)
    Btau_interval=int((Btau_num-1)/(Btau_show-1))
    Btau_choose=Btau_list[::Btau_interval]
    dirs = os.listdir('.')
    if not os.path.exists('./rhotheta'):
        os.mkdir('./rhotheta')
    files= os.listdir('./rhotheta')
    for file in files:
        if file.endswith('K_Btau.dat'):
            file_path = os.path.join('./rhotheta', file)
            os.remove(file_path)
    tmp_list = []
    for dir_name in dirs:
        if 'K.dat' in dir_name:
            number=re.findall(r'\d+\.\d+',dir_name)
            tmp_list.append(number[0])
    for temperature_name in tmp_list:
        temperature_name = str(temperature_name)
        extract_rho_theta('./{}K.dat'.format(temperature_name),dtheta=theta_interval,num_show=Btau_show, interval_line=Btau_interval,choose=Btau_choose,output_file='./{}K_Btau.dat'.format(temperature_name), output_dir='rhotheta')

