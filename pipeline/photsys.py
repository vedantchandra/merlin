def photsys():
        WAVE = {}
        WAVE['GALEX_FUV'] = [1539,255,515.0,'>','C4']
        WAVE['GALEX_NUV'] = [2316,729.9,781.8,'>','C4']

        WAVE['Tycho_B'] = [4208,686,4038.2,'s','C1']
        WAVE['Tycho_V'] = [5268,1037,3846.7,'s','C1']
        WAVE['Hipparcos_Hp'] = [5565.7,2388.5,3435.6,'s','C1']

        # WAVE['GaiaDR2_G']  = [6773.7,4358.4,2835.1,'+']
        # WAVE['GaiaDR2_BP'] = [5278.6,2279.4,3393.3,'+']
        # WAVE['GaiaDR2_RP'] = [7919.1,2943.7,2485.1,'+']

        WAVE['GaiaDR2_G']  = [6246.77,4358.43,2835.09,'+','C2']
        WAVE['GaiaDR2_BP'] = [5048.62,2279.45,3393.28,'+','C2']
        WAVE['GaiaDR2_RP'] = [7740.87,2943.72,2485.08,'+','C2']

        WAVE['GaiaEDR3_G']  = [6217.59,4052.97,3228.75,'+','C2']
        WAVE['GaiaEDR3_BP'] = [5109.71,2157.50,3552.01,'+','C2']
        WAVE['GaiaEDR3_RP'] = [7769.02,2924.44,2554.95,'+','C2']

        WAVE['GaiaMAW_G']   = [6229.36,4135.79,2847.56,'+','C2']
        WAVE['GaiaMAW_BPb'] = [5117.15,2434.57,3453.77,'+','C2']
        WAVE['GaiaMAW_BPf'] = [5037.40,2600.05,3412.20,'+','C2']
        WAVE['GaiaMAW_RP']  = [7752.07,3007.38,2482.27,'+','C2']

        WAVE['Bessell_U'] = [3605.1, 640.4,1803.1,'*','C7']
        WAVE['Bessell_B'] = [4400,   900.0,4000.0,'*','C7']
        WAVE['Bessell_V'] = [5512.1, 893.1,3579.8,'*','C7']
        WAVE['Bessell_R'] = [6585.9,1591.0,2971.4,'*','C7']
        WAVE['Bessell_I'] = [8059.9,1495.1,2405.3,'*','C7']

        # WAVE['SDSS_u'] = [3561.8,558.4,1568.5,'o']
        # WAVE['SDSS_g'] = [4718.9,1158.4,3965.9,'o']
        # WAVE['SDSS_r'] = [6185.2,1111.2,3162.0,'o']
        # WAVE['SDSS_i'] = [7499.7,1044.6,2602.0,'o']
        # WAVE['SDSS_z'] = [8961.5,1124.6,2244.7,'o']

        WAVE['SDSS_u'] = [3561.8,558.4 ,3631.0,'o','C3']
        WAVE['SDSS_g'] = [4718.9,1158.4,3631.0,'o','C3']
        WAVE['SDSS_r'] = [6185.2,1111.2,3631.0,'o','C3']
        WAVE['SDSS_i'] = [7499.7,1044.6,3631.0,'o','C3']
        WAVE['SDSS_z'] = [8961.5,1124.6,3631.0,'o','C3']

        WAVE['2MASS_J'] = [12350.0,1624.1,1594.0,'d','C5']
        WAVE['2MASS_H'] = [16620.0,2509.4,1024.0,'d','C5']
        WAVE['2MASS_Ks'] = [21590.0,2618.9,666.8,'d','C5']

        WAVE['WISE_W1'] = [33526.0 ,6626.4,309.5,'h','C6']
        WAVE['WISE_W2'] = [46028.0,10422.7,171.8,'h','C6']
        WAVE['WISE_W3'] = [115608.0,55055.7,31.7,'h','C6']
        WAVE['WISE_W4'] = [220883.0, 41016.8,8.4,'h','C6']

        # WAVE['PS_g'] = [4775.6,1166.5,3909.1,'x']
        # WAVE['PS_r'] = [6129.5,1318.1,3151.4,'x']
        # WAVE['PS_i'] = [7484.6,1242.6,2584.6,'x']
        # WAVE['PS_z'] = [8657.8,965.8,2273.1,'x']
        # WAVE['PS_y'] = [9603.1,614.9,2206.0,'x']

        WAVE['PS_g'] = [4775.6,1166.5,3631.0,'x','C0']
        WAVE['PS_r'] = [6129.5,1318.1,3631.0,'x','C0']
        WAVE['PS_i'] = [7484.6,1242.6,3631.0,'x','C0']
        WAVE['PS_z'] = [8657.8, 965.8,3631.0,'x','C0']
        WAVE['PS_y'] = [9603.1, 614.9,3631.0,'x','C0']

        return WAVE
def filtercurves():
        import numpy as np
        import socket
        hostname = socket.gethostname()
        if hostname[:4] == 'holy':
                RUNLOC = 'ODY'
        else:
                RUNLOC = 'LOCAL'

        if RUNLOC == 'LOCAL':
                filedir = '/Users/pcargile/Astro/gitrepos/h3/h3py/ms/filtercurves/'
        else:
                filedir = '/n/holyscratch01/conroy_lab/pacargile/ThePayne/filtercurves/'

        FC = {}
        for ff in photsys().keys():
                try:
                        FC[ff] = np.loadtxt(filedir+'{}.fc'.format(ff),dtype=[('wave',float),('trans',float)])
                except:
                        FC[ff] = {'wave':[],'trans':[]}
        return FC