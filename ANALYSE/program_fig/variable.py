import glob
type='t5'

if type=='t15':
    end_time=45
    active_time=16
    activation_time=30
    files=glob.glob('../donnee/t15_GMR_72F11.txt')
if type=='t5' :
    end_time=83
    active_time=39
    activation_time=45
    files=glob.glob('../donnee/t5_FCF_attP2.txt')

processus_nb=100
larva=1000