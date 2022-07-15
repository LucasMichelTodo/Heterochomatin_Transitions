#### Peak-Calling ####

os.makedirs('./Peak_Calling_MACS2', exist_ok=True)
outdir = './Peak_Calling_MACS2/'

params_form = ("-f BAMPE -B "
               "-g 2.41e7 "
               "--keep-dup all "
               "--fe-cutoff 1.5 "
               "--nomodel "
               "--extsize 150 "
               "-n {} "
               f"--outdir {outdir}")

for ip in IPs:
    prefix = ip.split('_')[0]
    inpt = [f for f in inputs if prefix in f][0]
    print(ip, inpt)
    pair = [ip, inpt]

    t = bamdir + pair[0]
    c = bamdir + pair[1]
    name = pair[0].split("_")[0]+'_'+pair[0].split("_")[1]+"_Macspeaks"
    params = params_form .format(name)

    macs2callpeak(t, c, params)

    print("==============================")
    print("Finished {}!" .format(name))
    print("==============================\n\n\n")
