#pragma once 
#include <vector>
#include <cmath> 
#include <tgmath.h>
#include <functional>
#include <iostream>
#include <array>
#include <algorithm>

#include "definitions.h"
#include "ranges.h"

#include <petscsys.h>


namespace numerical {

using namespace type;

struct operators {

    using ℤ = std::size_t;
    using ℝ = type::real_t; 

    operators(){};
    
    std::vector<ℝ> H(
        ℤ nodes, ℤ order=2, ℝ left=-1., ℝ right=1.);

    std::vector<ℝ> H_inverse(
        ℤ nodes, ℤ order=2, ℝ left=-1., ℝ right=1.);
};
    

/*
template<typename container> 
auto hinv_transform(container tail, std::size_t stop) -> long double(*)(long double){
    return [tail, stop](long double x){
        x = std::clamp(std::abs(x - (stop / 2)) - (stop / 2) + tail.size(), 0, tail.size());
        return stop * tail[x];
    };
}

template<std::size_t order> auto H_inverse(std::size_t N) {
    // std::array<long double, 1> xc    = {-1. , 1.}; 
    
    // std::array<long double, 3> d     = {-1/2, 0., 1/2}; 
    // std::array<long double, 3> bd    = {-1. , 1.}; 


    auto H_I = numerics::linrange<long double>(0., N, N + 1);

    H_I |= [N](long double x){

        constexpr std::array<long double, 1> tail = {1., 2.  }; 

        x = std::clamp(std::abs(x - (stop / 2)) - (stop / 2) + tail.size(), 0, tail.size());
        return stop * tail[x];
    };
    
}
*/
/* 
function D1(p, N; xc = (-1, 1))

  if p == 2
    bhinv = [2]
    d  = [-1/2 0 1/2]
    bd  = [-1 1]
  elseif p == 4
    bhinv = [48/17 48/59 48/43 48/49]
    d  = [1/12 -2/3 0 2/3 -1/12]
    bd = [-24/17  59/34  -4/17  -3/34  0     0;
           -1/2    0      1/2    0     0     0;
            4/43 -59/86   0     59/86 -4/43  0;
            3/98   0    -59/98   0    32/49 -4/49]
  elseif p == 6
    bhinv = [43200/13649 8640/12013 4320/2711 4320/5359 8640/7877 43200/43801];

    d = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60];

    x1=0.70127127127127;

    bd = [-21600/13649              8(16200x1-953)/40947     (-1036800x1+715489)/81894 3(86400x1-62639)/13649    5(-207360x1+147127)/81894 (129600x1-89387)/40947     0           0           0;
          8(-16200x1+953)/180195    0                        (86400x1-57139)/12013     (-1036800x1+745733)/72078 5(25920x1-18343)/12013    (-345600x1+240569)/120130  0           0           0;
          (1036800x1-715489)/162660 (-86400x1+57139)/5422    0                         (259200x1-176839)/8133    (-345600x1+242111)/10844  (259200x1-182261)/27110    0           0           0;
          3(-86400x1+62639)/53590   (1036800x1-745733)/64308 (-259200x1+176839)/16077  0                         (259200x1-165041)/32154   (-1036800x1+710473)/321540 72/5359     0           0;
          (207360x1-147127)/47262   5(-25920x1+18343)/7877   (345600x1-242111)/15754   (-259200x1+165041)/23631  0                         8640x1/7877                -1296/7877  144/7877    0;
          (-129600x1+89387)/131403  (345600x1-240569)/87602  (-259200x1+182261)/43801  (1036800x1-710473)/262806 -43200x1/43801            0                          32400/43801 -6480/43801 720/43801];

  elseif p == 8
    bhinv = [5080320/1498139 725760/1107307 80640/20761 725760/1304999 725760/299527 80640/103097 725760/670091 5080320/5127739];
    d = [1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280];

    bd = [    -2540160/1498139        699846290793/311403172540  -10531586157/311403172540  -145951651079/186841903524    398124597/15570158627    39152001/113858564  -80631892229/934209517620  -6230212503/311403172540        0                 0               0               0;
          -24132630717/55557028660               0                 2113176981/23016483302      5686186719/11508241651   -3408473341/138098899812  -39291999/210388330     607046586/11508241651    3460467023/483346149342        0                 0               0               0;
            3510528719/90623010660      -704392327/1294614438               0                   503511235/1294614438      354781619/2589228876       407439/3944590        -2986842/16597621        169381493/3020767022          0                 0               0               0;
          145951651079/1139279786988   -5686186719/13562854607    -1510533705/27125709214               0                6763379967/54251418428    13948923/49589962    -1603900430/40688563821   -3742312557/189879964498        0                 0               0               0;
            -398124597/21790888777      3408473341/37355809332    -1064344857/12451936444     -6763379967/12451936444             0                  763665/1198108     -1282435899/12451936444    7822226819/261490665324    -2592/299527          0               0               0;
              -1864381/23506116           13097333/58765290           -407439/19588430           -4649641/11753058          -254555/1237164               0                 5346432/9794215           -923328/9794215          3072/103097       -288/103097        0               0;
           11518841747/417855345780     -607046586/6964255763         1192698/23768791         1603900430/20892767289    1282435899/27857023052   -48117888/63658645              0                 301190400/366539777     -145152/670091      27648/670091    -2592/670091        0;
            6230212503/1065851828540   -3460467023/319755548562   -1524433437/106585182854     3742312557/106585182854  -7822226819/639511097124   58169664/487135205   -2108332800/2804873233              0               4064256/5127739  -1016064/5127739  193536/5127739  -18144/5127739];
  elseif p == 10
    bhinv = [18289152000/5261271563 1828915200/2881040311 406425600/52175551 6096384/11662993 87091200/50124587 72576000/50124587 87091200/148333439 152409600/63867949 16257024/20608675 1828915200/1704508063 18289152000/18425967263];

    d = [-1/1260,5/504,-5/84,5/21,-5/6,0,5/6,-5/21,5/84,-5/504,1/1260];

    bd = [  -1.7380923775745425e+00   2.3557601935237220e+00  -1.5328406598563976e-01  -5.7266565770416333e-01  -1.8308103515008173e-01   1.8186748267946842e-01   2.0034232582598244e-01   2.2678007363666621e-02  -1.1782459320459637e-01  -3.0591175636402144e-02   3.4890895862586133e-02   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
            -4.3020203737210871e-01   0.0000000000000000e+00   1.1837297346927406e-01   3.3928601158526644e-01   1.3241927733034406e-01  -8.7495003780608913e-02  -1.1750484124279399e-01  -1.6401912273575153e-02   6.2537843443041474e-02   1.7143274696828435e-02  -1.8155585855667674e-02   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
             3.4348531361887280e-01  -1.4525207124434036e+00   0.0000000000000000e+00   2.9011513992277767e+00  -2.2419288742360557e+00  -5.4662873578741478e-01   1.2908050607446131e+00   6.1514504292452719e-02  -4.2442625460011202e-01   1.5579158905288801e-02   5.2969140277981920e-02   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
             8.6111387816878188e-02  -2.7937273515056432e-01  -1.9467880944770807e-01   0.0000000000000000e+00   2.0170150914578375e-01   2.4269917331475005e-01  -7.7261988327590472e-02   5.0649247607525059e-02  -7.4775049946661561e-03  -4.0978487203372188e-02   1.8608207238964152e-02   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
             9.1509035082611684e-02  -3.6243526359648576e-01   5.0007055839856984e-01  -6.7045605191055857e-01   0.0000000000000000e+00  -1.7807807859119628e-02   7.5000761407401195e-01  -2.2979723229714316e-01  -1.2521154324370892e-01   6.8278284106004450e-02  -4.1575927541817690e-03   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
            -7.5752056274147259e-02   1.9956355926115746e-01   1.0160630736447970e-01  -6.7227694623145351e-01   1.4839839882599690e-02   0.0000000000000000e+00   5.4091068834671807e-01  -1.2712520372174399e-01  -8.9292453564020990e-02   1.6181541970619609e-01  -5.4289154769785249e-02   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
            -3.3838029883391296e-02   1.0867927550524317e-01  -9.7293058702223670e-02   8.6783825404790446e-02  -2.5344131542932297e-01  -2.1934035945002228e-01   0.0000000000000000e+00   2.7184438867288430e-01   1.9102691945078512e-01  -4.8646826827046824e-02  -6.2407959378425991e-03   4.6597719614658163e-04   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
            -1.5567948806367624e-02   6.1656604470023607e-02  -1.8844858059892756e-02  -2.3122780265804038e-01   3.1560994521078772e-01   2.0951677187991255e-01  -1.1048784865195491e+00   0.0000000000000000e+00   1.1823059621092409e+00  -5.3610400867086083e-01   1.5931375952374752e-01  -2.3673846172827626e-02   1.8939076938262100e-03   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00;
             2.6737701764454301e-02  -7.7712278574126673e-02   4.2981266272823705e-02   1.1284579710276557e-02   5.6847566375570611e-02   4.8647834370398067e-02  -2.5665536068472994e-01  -3.9083324869946684e-01   0.0000000000000000e+00   6.5716944195909766e-01  -1.5822272208022428e-01   4.6954983762905661e-02  -7.8258306271509429e-03   6.2606645017207550e-04   0.0000000000000000e+00   0.0000000000000000e+00;
             9.4425181052687698e-03  -2.8976375375532045e-02  -2.1459742428921558e-03   8.4117843695442701e-02  -4.2165149106440383e-02  -1.1991463562335723e-01   8.8902467992349743e-02   2.4105392677971343e-01  -8.9388344421253152e-01   0.0000000000000000e+00   8.6496680152924643e-01  -2.5547312415382800e-01   6.3868281038457000e-02  -1.0644713506409501e-02   8.5157708051276015e-04   0.0000000000000000e+00;
            -9.9625965676187218e-03   2.8387641187789508e-02  -6.7495090936003027e-03  -3.5335033597892078e-02   2.3750992019053968e-03   3.7216380474824604e-02   1.0550378667904333e-02  -6.6265458456725809e-02   1.9908619649258188e-01  -8.0014409359906680e-01   0.0000000000000000e+00   8.2714572225493910e-01  -2.3632734921569687e-01   5.9081837303924217e-02  -9.8469728839873684e-03   7.8775783071898962e-04];
  else
    error(string("Operators for order ", p, " are not implemented"))
  end

  (bm, bn) = size(bd);
  Np = N+1

  if(Np < 2*bm || Np < bn)
    error("Grid not big enough to support the operator. Grid must have N >= ", max(bn,2*bm))
  end

  h = (xc[2] - xc[1]) / N
  @assert h > 0
  H_I = 1:Np
  H_V = ones(Np)
  H_V[1:bm] = 1 ./ bhinv[:]
  H_V[Np-bm+1:Np] = 1 ./ bhinv[end:-1:1]
  H = sparse(H_I, H_I, h * H_V)
  HI = sparse(H_I, H_I, 1 ./ (h * H_V))

  n = floor(Int64, length(d)/2);
  B_I1 = (bm+1:Np-bm) * ones(1,p+1)
  B_J1 = ones(Np-2bm,1) * (-div(p,2):div(p,2))' + B_I1
  B_V1 = ones(Np-2bm)
  B_V1 = kron(B_V1, d)

  B_I2 = (1:bm) * ones(1, bn)
  B_J2 = ones(bm) * (1:bn)'
  B_V2 = bd
  B_I3 = (Np+1) .- B_I2
  B_J3 = (Np+1) .- B_J2
  B_V3 = -B_V2
  D = sparse([B_I1[:];B_I2[:];B_I3[:]],
             [B_J1[:];B_J2[:];B_J3[:]],
             [B_V1[:];B_V2[:];B_V3[:]]/h, N+1, N+1)

  r = range(xc[1], stop=xc[2], length=N+1)

  (D, HI, H, r)
end

*/

}; /* numerical:: */