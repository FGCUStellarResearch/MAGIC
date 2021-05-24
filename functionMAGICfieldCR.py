from __future__ import print_function, division
from tabulate import tabulate
def MAGIC_Field_CR(ra_pos, dec_pos):
    import numpy as np 
    from astroquery.vizier import Vizier
    import astropy.units as u
    import astropy.coordinates as coord
    from PyAstronomy import pyasl
    import matplotlib.pyplot as plt
    #from tabulate import tabulate
    from astropy.table import Table, Column, MaskedColumn, QTable
    from quantiphy import Quantity
    from astropy.io import ascii
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(coord.SkyCoord(ra=ra_pos, dec=dec_pos,unit=(u.deg, u.deg),frame='icrs'),width="7d",catalog="V/125")
    if result ==[]:
        data=np.array([ra_pos, dec_pos, 0, 0, 0])
        #data = data.transpose
        return data
    t = result['V/125/obcat']
    #t.columns
    #use ALS IDs to ensure uniqueness


    ALS0 = np.array(t['ALS'])
    ra = np.array(t['RAJ2000'])
    dec = np.array(t['DEJ2000'])
    #ALS0[0:3]
    t.columns  
    
    #Section 3 of Jupyter Notebook 

    from astroquery.simbad import Simbad
    #Right now that's B, V, sptype, but later might want more...
    customSimbad = Simbad()
    customSimbad.add_votable_fields('sptype','fe_h','flux(V)','flux(B)','flux(U)')
    vmag = np.zeros(len(ALS0))
    bmag = np.zeros(len(ALS0))
    st = np.empty(len(ALS0), dtype='S2')
    #vmag = np.zeros(200)
    #bmag = np.zeros(200)
    for ii in range(len(ALS0)):
        result_table = customSimbad.query_object('ALS'+np.array2string(ALS0[ii]))
        #print(str(result_table[0][11])[0:2])
        #print(ALS0[ii])
        
        
        if result_table is None:
            print('None 2')
            vmag[ii] = np.nan
            bmag[ii] = np.nan
            st[ii] = 'XX'
        else:
            vmag[ii] = np.asarray(result_table['FLUX_V'])
            bmag[ii] = np.asarray(result_table['FLUX_B'])
            st[ii] = str(result_table[0][11])[0:2]

    print(bmag)
    idx = (np.isnan(vmag))
    print(idx)

    vmag = vmag[~idx]
    bmag = bmag[~idx]
    st = st[~idx]
    ra = ra[~idx]
    dec = dec[~idx]

    idx = (np.isnan(bmag))
    print(idx)

    vmag = vmag[~idx]
    bmag = bmag[~idx]
    st = st[~idx]
    ra = ra[~idx]
    dec = dec[~idx]

    print(bmag)  

    #Start of section 4 in Jupyter notebook 
    #collect spectral type, V magnitude information for each star

    #Use SIMBAD spectral types, recorded in st string array
    #Do some temporary cleanup here. For stars with generic OB type, assign type O5V. For Wolf-Rayet stars (W..), 
    #also assign type O5V. Note that this is definitely too cool, but right now don't have the models to assign a hotter type
      
    st_vals = list()
    for ii in range(len(st)):
        temp1 = np.array2string(st[ii])[2:4]
        if len(temp1)== 2 and temp1[1]=="'":
         temp1 = temp1[0]+'5'
        if temp1 =="OB":
            temp1 = "O5"
        if temp1[0] == "W":
            temp1 = "O5"
        if temp1[0]=="'":
            temp1="O5"
        if len(temp1)== 1 and temp1[1]=="'":
            temp1="G5"

        temp1 = temp1+'V'
        st_vals.append(temp1)
    st_vals = np.asarray(st_vals)
    len(st_vals)
    #MeanStars is a python package which includes Mamajek mean stars catalog
    from MeanStars import MeanStars
    ms = MeanStars()
    ms_st = ms.data

    #ms_st are spectral types from Mamajek mean stars catalog
    ms_st = ms_st['SpT']
    #ms_temp = ms_st['Teff']
    #ms_temp are temperatures from Mamajek catalog
    ms_temp = ms.data['Teff']
    #match up rows in Mamajek catalog with observed STs
    sloc = []
    for x in range(len(st_vals)):
        sloc_temp = np.nonzero(st_vals[x]==ms_st)
        sloc.append(sloc_temp[0])

    #star_temp = temperature of stars suitable for calling to function
    #Note for stars with no temperature value, assign T<10,000K so they will be ignored later
    #For stars with T>40000, assign 40000 (model fits for >40,000 are broken for now)

    star_temp = []
    for x in range(len(st_vals)):
        star_temp.append(ms_temp[sloc[x]])
        if len(star_temp[x]) == 0:
            star_temp[x]=9999
        if star_temp[x]>40000:
            star_temp[x]=40000

    #Estimate count rate for each star based on spectral type, magnitude, model atmosphere fits
    #Polynomial fits to the atmosphere fits derived from my MATLAB code (which needs to be ported to Python!)
    #cr1   -0.3375  -34.1223  227.7543 -172.6182   34.0610
    #cr2 3.2407  -51.6510  219.8193  -91.3053    5.2321

    cr1 = np.zeros(len(st_vals))
    cr2 = np.zeros(len(st_vals))
    for x in range(len(st_vals)):
        cr2[x] = (2.512**-bmag[x])*1e4*np.polyval([-0.3375,-34.1223,227.7543,-172.6182,34.0610],1e-4*star_temp[x])    
        cr1[x] = (2.512**-bmag[x])*1e4*np.polyval([3.2407,-51.6510,219.8193,-91.3053,5.2321],1e-4*star_temp[x])
        
    #plt.figure(1)
    #plt.plot(star_temp,cr1,'.')
    #plt.show()      
    #Summed total count rates over field in both bands, as well as max local CRs in field in both bands
    print(cr1)
    tot_cr1 = np.nansum(cr1)
    tot_cr2 = np.nansum(cr2)
    
    count1 = np.count_nonzero(cr1 > 100)
    count2 = np.count_nonzero(cr2 > 100)
    #np.nanmax(cr1)
    #np.nanmax(cr2)

    np.nanmedian(cr1)
    #np.nanmedian(cr2)

    #Histogram of count rates in field
    #plt.figure(2)
    #plt.hist(np.log10(cr1[~np.isnan(cr1)]),bins=40)
    #plt.show()


    #Scatter plot
    #pyastro and future import was moved to the beginning of the file
    #plt.plot(ra,dec,'x')
    #plt.show()


    ra_d = np.zeros(len(ra))
    dec_d = np.zeros(len(dec))
    for ii in range(len(ra)):
        temp = str(ra[ii])
        temp2 = str(dec[ii])
        ra_d[ii], dec_d[ii] = pyasl.coordsSexaToDeg(temp[0:-1]+" "+temp2[0:-1])
    
    #plt.figure(3)
    #plt.scatter(ra_d,dec_d,cr1*0.01,alpha=0.3)
    #plt.axis('square')
    #plt.xlabel('RA')
    #plt.ylabel('Dec')
    #plt.show()

    #Number of stars with T>10000
    tstar = np.asarray(star_temp)
    len(tstar[tstar>10000]) 
    
     
    b_i, v_i = cr1, cr2
    
    #print(b_i, v_i) 
    maxcr1= np.max(b_i)
    #x=f"{maxcr1.value:0.03f}"
    maxcr2= np.max(v_i)
    in_cr=maxcr1+maxcr2
    data=np.array([ra_pos, dec_pos, maxcr1,maxcr2,tot_cr1, tot_cr2, count1,count2])
    #data = data.transpose
    #print(data)
    #print(maxcr2)
    #data=(np.array(x,dtype=object,copy=False, order='C'))
    #fopen = open("sample.txt", "a") 
    #print(fopen)
    #np.savetxt("sample.txt", data, fmt='%10.5f', newline=' ')
    
    #np.savetxt("Datamagic3",data, fmt='%10.5f',newline=' ')
   # data=(np.array(x,dtype=object,copy=False, order='C')
    #ascii.write(data, 'values.dat', names=['x'], overwrite=True)   
    #ascii.write(data,'values.csv',format='csv',fast_writer=False) 
   # print(x)
   # y = maxcr2
    
    #data1=QTable([maxcr1], names=['Maxcr1'])
    
    #lines = ['in_CR                   & maxCr1            & maxCr2       ','----------------------- & ----------------- & -------------','              32876.33701130775 & 18225.588784508154 & 14650.748916622626']
   # data = ascii.read(lines, data_start=2, delimiter='&')
   # print(data)
  # plt.hist(np.log10(cr1[~np.isnan(cr1)]),bins=40)
 #   return plt.hist(np.log10(cr1[~np.isnan(cr1)]),bins=40)
 
    return data
    


#MAGIC_Field_CR(300,33)

