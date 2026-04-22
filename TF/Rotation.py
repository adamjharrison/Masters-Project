import numpy as np

rotmat = np.array([[1.0,0.0,0.0],[0.0,0.0,-1.0],[0.0,1.0,0.0]]) #Rotates 90deg wrt x axis

tf = np.array([[-8.548130057571822E-03,2.719924568325373E-03,1.727290824601101E+00,-1.711796134105602E+00],
               [-1.343054908593903E-02,1.304004028490127E-03,5.927355424533597E-02,5.926045976573367E-02],
               [7.359856293247700E-02,3.042380162073215E+00,-9.851492698813099E-01,-1.128779753990537E+00]])
tf_mom = np.array([[-3.644383604201998E+00,9.317393615136383E-01,4.708983449820484E-01,2.241745897706311E+00],
                  [4.107341723179617E+00,-1.065102626702020E+00,-1.521128122458516E+00,-1.521110974019082E+00],
                  [3.615046233265347E-02,1.749642388052235E+00,-9.088314593095068E-01,-8.769613910753812E-01]])

#Applies rotation matrix onto position and momentum coordinates
tf_flipx = rotmat@tf
tf_mom_flipx = rotmat@tf_mom

np.set_printoptions(precision=17, suppress=False)
print("Flipped TF Structure:\n")
print(tf_flipx)
print("Flipped Momentum of TF Structure:\n")
print(tf_mom_flipx)