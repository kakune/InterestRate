/**
 * @file integral_1d_adaptive.tpp
 * @brief This implements the adaptive 1d-integral functions
 * @author kakune
 * @date 5/3/2024
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <queue>
#include <tuple>
#include <vector>
#ifdef NINCLUDE_TPP
#include "math/integral_1d.hpp"
#endif

namespace Math::Integral
{
namespace FiniteInterval
{

#define MATH_INTEGRAL_ADAPTIVE_DEFAULT_MAX_SIZE_INTERVAL 1000
#define MATH_INTEGRAL_ADAPTIVE_DEFAULT_MAX_EVAL_FUNC 10000
constexpr double gCoeffC    = 32.0;
constexpr double gCoeffD    = 32.0;
constexpr double gFactorEps = 50.0 * std::numeric_limits<double>::epsilon();

template <std::size_t NPoints_, typename TypeY_>
std::tuple<TypeY_, double, double, double> evalNewtonCotesIntegralAndErr(
    double inLeft, double inRight, const std::vector<TypeY_>& inVals );

template <std::size_t MaxNInterval_ =
              MATH_INTEGRAL_ADAPTIVE_DEFAULT_MAX_SIZE_INTERVAL>
class Errors;

template <typename Func_, typename TypeY_ = double,
          std::size_t MaxNInterval_ =
              MATH_INTEGRAL_ADAPTIVE_DEFAULT_MAX_SIZE_INTERVAL,
          std::size_t MaxNFuncEval_ =
              MATH_INTEGRAL_ADAPTIVE_DEFAULT_MAX_EVAL_FUNC>
class DoublyAdaptiveNewtonCotes;

constexpr std::array<double, 5> gWeightForCotes5 = {
    1.5555555555555556e-01, 7.1111111111111114e-01, 2.6666666666666666e-01,
    7.1111111111111114e-01, 1.5555555555555556e-01,
};
constexpr std::array<double, 20> gNullWeightForCotes5 = {
    1.2710311885185779e-01,  -5.0841247540743117e-01, 7.6261871311114682e-01,
    -5.0841247540743117e-01, 1.2710311885185779e-01,  -3.3628324334270127e-01,
    6.7256648668540253e-01,  0.0000000000000000e+00,  -6.7256648668540253e-01,
    3.3628324334270127e-01,  5.6842242780997809e-01,  -2.8421121390498905e-01,
    -5.6842242780997809e-01, -2.8421121390498905e-01, 5.6842242780997809e-01,
    -6.7256648668540253e-01, -3.3628324334270127e-01, 0.0000000000000000e+00,
    3.3628324334270127e-01,  6.7256648668540253e-01,
};
constexpr std::array<double, 9> gWeightForCotes9 = {
    6.9770723104056437e-02,  4.1537918871252205e-01,  -6.5467372134038804e-02,
    7.4045855379188708e-01,  -3.2028218694885363e-01, 7.4045855379188708e-01,
    -6.5467372134038804e-02, 4.1537918871252205e-01,  6.9770723104056437e-02,
};
constexpr std::array<double, 72> gNullWeightForCotes9 = {
    1.1018547692345270e-02,  -8.8148381538762158e-02, 3.0851933538566756e-01,
    -6.1703867077133512e-01, 7.7129833846416884e-01,  -6.1703867077133512e-01,
    3.0851933538566756e-01,  -8.8148381538762158e-02, 1.1018547692345270e-02,
    -4.2674651711845396e-02, 2.5604791027107238e-01,  -5.9744512396583560e-01,
    5.9744512396583560e-01,  0.0000000000000000e+00,  -5.9744512396583560e-01,
    5.9744512396583560e-01,  -2.5604791027107238e-01, 4.2674651711845396e-02,
    1.1236757938944257e-01,  -4.7756221240513091e-01, 6.1802168664193424e-01,
    2.8091894847360643e-02,  -5.6183789694721287e-01, 2.8091894847360643e-02,
    6.1802168664193424e-01,  -4.7756221240513091e-01, 1.1236757938944257e-01,
    -2.3112700627433047e-01, 6.3559926725440874e-01,  -2.3112700627433047e-01,
    -5.2003576411724362e-01, 0.0000000000000000e+00,  5.2003576411724362e-01,
    2.3112700627433047e-01,  -6.3559926725440874e-01, 2.3112700627433047e-01,
    3.9111964345081668e-01,  -5.8667946517622505e-01, -3.0730829128278453e-01,
    2.5143405650409645e-01,  5.0286811300819290e-01,  2.5143405650409645e-01,
    -3.0730829128278453e-01, -5.8667946517622505e-01, 3.9111964345081668e-01,
    -5.5619114160254801e-01, 2.7809557080127401e-01,  5.1646320291665171e-01,
    3.5755144817306656e-01,  0.0000000000000000e+00,  -3.5755144817306656e-01,
    -5.1646320291665171e-01, -2.7809557080127401e-01, 5.5619114160254801e-01,
    6.6477556470172228e-01,  1.6619389117543057e-01,  -1.8993587562906350e-01,
    -4.0361373571175990e-01, -4.7483968907265878e-01, -4.0361373571175990e-01,
    -1.8993587562906350e-01, 1.6619389117543057e-01,  6.6477556470172228e-01,
    -6.4550259924248821e-01, -4.8412694943186613e-01, -3.2275129962124410e-01,
    -1.6137564981062205e-01, 0.0000000000000000e+00,  1.6137564981062205e-01,
    3.2275129962124410e-01,  4.8412694943186613e-01,  6.4550259924248821e-01,
};
constexpr std::array<double, 17> gWeightForCotes17 = {
    3.0797894233299011e-02,  2.6128238288028033e-01,  -3.6795289329867603e-01,
    1.7037379778090087e+00,  -3.9501480717783930e+00, 8.5525299934402952e+00,
    -1.3934614237197881e+01, 1.9180342211078734e+01,  -2.0951950514333333e+01,
    1.9180342211078734e+01,  -1.3934614237197881e+01, 8.5525299934402952e+00,
    -3.9501480717783930e+00, 1.7037379778090087e+00,  -3.6795289329867603e-01,
    2.6128238288028033e-01,  3.0797894233299011e-02,
};
constexpr std::array<double, 255> gNullWeightForCotes17 = {
    1.7047365484036969e-03,  -2.7275784774459150e-02, 2.0456838580844364e-01,
    -9.5465246710607032e-01, 3.1026205180947288e+00,  -7.4462892434273487e+00,
    1.3651530279616807e+01,  -1.9502186113738293e+01, 2.1939959377955578e+01,
    -1.9502186113738293e+01, 1.3651530279616807e+01,  -7.4462892434273487e+00,
    3.1026205180947288e+00,  -9.5465246710607032e-01, 2.0456838580844364e-01,
    -2.7275784774459150e-02, 1.7047365484036969e-03,  -9.4915714022159627e-03,
    1.3288199963102348e-01,  -8.5424142619943655e-01, 3.3220499907755867e+00,
    -8.6373299760165256e+00, 1.5547193956829744e+01,  -1.9002125947236355e+01,
    1.3572947105168828e+01,  0.0000000000000000e+00,  -1.3572947105168828e+01,
    1.9002125947236355e+01,  -1.5547193956829744e+01, 8.6373299760165256e+00,
    -3.3220499907755867e+00, 8.5424142619943655e-01,  -1.3288199963102348e-01,
    9.4915714022159627e-03,  3.6721149063598142e-02,  -4.4524393239612747e-01,
    2.4144155509315781e+00,  -7.6150482870636651e+00, 1.5037310541543439e+01,
    -1.7961232035732444e+01, 9.1894675531654357e+00,  7.2202959346299851e+00,
    -1.5753372948283605e+01, 7.2202959346299851e+00,  9.1894675531654357e+00,
    -1.7961232035732444e+01, 1.5037310541543439e+01,  -7.6150482870636651e+00,
    2.4144155509315781e+00,  -4.4524393239612747e-01, 3.6721149063598142e-02,
    -1.1389885682659154e-01, 1.1817006395758871e+00,  -5.2962968424365071e+00,
    1.3027181749541409e+01,  -1.7768221664948282e+01, 9.4393677595037744e+00,
    8.1437682631012951e+00,  -1.4251594460427265e+01, 0.0000000000000000e+00,
    1.4251594460427265e+01,  -8.1437682631012951e+00, -9.4393677595037744e+00,
    1.7768221664948282e+01,  -1.3027181749541409e+01, 5.2962968424365071e+00,
    -1.1817006395758871e+00, 1.1389885682659154e-01,  3.0014984205120593e-01,
    -2.6263111179480521e+00, 9.4697275167155475e+00,  -1.7063518520611058e+01,
    1.2696338318766013e+01,  5.4777346174345078e+00,  -1.4032005115893879e+01,
    -1.1555768918971430e+00, 1.3866922702765715e+01,  -1.1555768918971430e+00,
    -1.4032005115893879e+01, 5.4777346174345078e+00,  1.2696338318766013e+01,
    -1.7063518520611058e+01, 9.4697275167155475e+00,  -2.6263111179480521e+00,
    3.0014984205120593e-01,  -6.9333963835801982e-01, 5.0267123780956435e+00,
    -1.4040127676749902e+01, 1.6466816411002974e+01,  -6.9333963835801982e-01,
    -1.4040127676749902e+01, 1.9066840054845546e+00,  1.3346788038391882e+01,
    0.0000000000000000e+00,  -1.3346788038391882e+01, -1.9066840054845546e+00,
    1.4040127676749902e+01,  6.9333963835801982e-01,  -1.6466816411002974e+01,
    1.4040127676749902e+01,  -5.0267123780956435e+00, 6.9333963835801982e-01,
    1.4311829358292953e+00,  -8.4081997479971111e+00, 1.7174195229951543e+01,
    -9.5327006261487011e+00, -1.0938326723838186e+01, 7.8970629852009342e+00,
    1.1858372896871305e+01,  -3.0412637386372525e+00, -1.2880646422463659e+01,
    -3.0412637386372525e+00, 1.1858372896871305e+01,  7.8970629852009342e+00,
    -1.0938326723838186e+01, -9.5327006261487011e+00, 1.7174195229951543e+01,
    -8.4081997479971111e+00, 1.4311829358292953e+00,  -2.6735921944284664e+00,
    1.2365363899231658e+01,  -1.6709951215177917e+01, -1.6709951215177918e+00,
    1.3367960972142335e+01,  6.3497814617676083e+00,  -8.6891746318925165e+00,
    -1.1696965850624542e+01, 0.0000000000000000e+00,  1.1696965850624542e+01,
    8.6891746318925165e+00,  -6.3497814617676083e+00, -1.3367960972142335e+01,
    1.6709951215177918e+00,  1.6709951215177917e+01,  -1.2365363899231658e+01,
    2.6735921944284664e+00,  4.5591565530058782e+00,  -1.5957047935520574e+01,
    1.1397891382514695e+01,  1.1397891382514695e+01,  -4.3838043778902680e+00,
    -1.2800708783439580e+01, -6.4880304792775965e+00, 6.1373261290463752e+00,
    1.2274652258092750e+01,  6.1373261290463752e+00,  -6.4880304792775965e+00,
    -1.2800708783439580e+01, -4.3838043778902680e+00, 1.1397891382514695e+01,
    1.1397891382514695e+01,  -1.5957047935520574e+01, 4.5591565530058782e+00,
    -7.1376364508808008e+00, 1.7844091127202002e+01,  -2.1412909352642404e+00,
    -1.3561509256673522e+01, -8.1808294706249178e+00, 4.1178671832004623e+00,
    1.1804552591841325e+01,  9.6083567608010778e+00,  0.0000000000000000e+00,
    -9.6083567608010778e+00, -1.1804552591841325e+01, -4.1178671832004623e+00,
    8.1808294706249178e+00,  1.3561509256673522e+01,  2.1412909352642404e+00,
    -1.7844091127202002e+01, 7.1376364508808008e+00,  1.0294045683708612e+01,
    -1.6727824236026496e+01, -7.7205342627814595e+00, 6.4337785523178832e+00,
    1.2669594687641368e+01,  9.2052523902394316e+00,  1.9796241699439637e-01,
    -8.4134027222618464e+00, -1.1877745019663783e+01, -8.4134027222618464e+00,
    1.9796241699439637e-01,  9.2052523902394316e+00,  1.2669594687641368e+01,
    6.4337785523178832e+00,  -7.7205342627814595e+00, -1.6727824236026496e+01,
    1.0294045683708612e+01,  -1.3692360757232027e+01, 1.1980815662578024e+01,
    1.3692360757232027e+01,  5.1346352839620097e+00,  -4.7396633390418552e+00,
    -1.0927557142790944e+01, -1.1585843717657870e+01, -7.2411523235361672e+00,
    0.0000000000000000e+00,  7.2411523235361672e+00,  1.1585843717657870e+01,
    1.0927557142790944e+01,  4.7396633390418552e+00,  -5.1346352839620097e+00,
    -1.3692360757232027e+01, -1.1980815662578024e+01, 1.3692360757232027e+01,
    1.6769648614663378e+01,  -4.1924121536658445e+00, -1.2577236460997533e+01,
    -1.2577236460997533e+01, -7.7398378221523281e+00, -9.6747972776904101e-01,
    5.4823851240245656e+00,  9.9972905202800906e+00,  1.1609756733228490e+01,
    9.9972905202800906e+00,  5.4823851240245656e+00,  -9.6747972776904101e-01,
    -7.7398378221523281e+00, -1.2577236460997533e+01, -1.2577236460997533e+01,
    -4.1924121536658445e+00, 1.6769648614663378e+01,  -1.8797050108382589e+01,
    -4.6992625270956472e+00, 4.6992625270956472e+00,  1.0069848272347816e+01,
    1.2083817926817378e+01,  1.1412494708660859e+01,  8.7272018360347730e+00,
    4.6992625270956472e+00,  0.0000000000000000e+00,  -4.6992625270956472e+00,
    -8.7272018360347730e+00, -1.1412494708660859e+01, -1.2083817926817378e+01,
    -1.0069848272347816e+01, -4.6992625270956472e+00, 4.6992625270956472e+00,
    1.8797050108382589e+01,  1.8987887997058074e+01,  1.1867429998161299e+01,
    5.6963663991174229e+00,  4.7469719992645187e-01,  -3.7975775994116150e+00,
    -7.1204579988967787e+00, -9.4939439985290370e+00, -1.0918035598308393e+01,
    -1.1392732798234846e+01, -1.0918035598308393e+01, -9.4939439985290370e+00,
    -7.1204579988967787e+00, -3.7975775994116150e+00, 4.7469719992645187e-01,
    5.6963663991174229e+00,  1.1867429998161299e+01,  1.8987887997058074e+01,
};
constexpr std::array<double, 33> gWeightForCotes33 = {
    1.4695535312430933e-02,  1.4401180552105430e-01,  -3.0306061552847280e-01,
    1.5082956203676328e+00,  -4.4611013899969842e+00, 1.1036738075229430e+01,
    -2.0246745867038545e+01, 2.7624501654321350e+01,  -2.4354757093571234e+01,
    7.4467704632794574e+00,  1.3501113227549911e+01,  -1.9331291092852421e+01,
    3.9483050394538197e+00,  1.5326402464090805e+01,  -1.4231233070221510e+01,
    -5.4208814966182697e+00, 1.7596473481403098e+01,  -5.4208814966182697e+00,
    -1.4231233070221510e+01, 1.5326402464090805e+01,  3.9483050394538197e+00,
    -1.9331291092852421e+01, 1.3501113227549911e+01,  7.4467704632794574e+00,
    -2.4354757093571234e+01, 2.7624501654321350e+01,  -2.0246745867038545e+01,
    1.1036738075229430e+01,  -4.4611013899969842e+00, 1.5082956203676328e+00,
    -3.0306061552847280e-01, 1.4401180552105430e-01,  1.4695535312430933e-02,
};
constexpr std::array<double, 495> gNullWeightForCotes33 = {
    5.8715078831162645e-08,  -1.8788825225972046e-06, 2.9122679100256674e-05,
    -2.9122679100256673e-04, 2.1113942347686090e-03,  -1.1823807714704208e-02,
    5.3207134716168940e-02,  -1.9762650037434179e-01, 6.1758281366981804e-01,
    -1.6468875031195149e+00, 3.7878412571748843e+00,  -7.5756825143497686e+00,
    1.3257444400112094e+01,  -2.0396068307864759e+01, 2.7680378417816460e+01,
    -3.3216454101379753e+01, 3.5292482482715990e+01,  -3.3216454101379753e+01,
    2.7680378417816460e+01,  -2.0396068307864759e+01, 1.3257444400112094e+01,
    -7.5756825143497686e+00, 3.7878412571748843e+00,  -1.6468875031195149e+00,
    6.1758281366981804e-01,  -1.9762650037434179e-01, 5.3207134716168940e-02,
    -1.1823807714704208e-02, 2.1113942347686090e-03,  -2.9122679100256673e-04,
    2.9122679100256674e-05,  -1.8788825225972046e-06, 5.8715078831162645e-08,
    4.6603649039042801e-07,  -1.3981094711712842e-05, 2.0225983682944576e-04,
    -1.8781270562734248e-03, 1.2569004145829843e-02,  -6.4520887948593189e-02,
    2.6394908706242670e-01,  -8.8234409103725508e-01, 2.4509558084368197e+00,
    -5.7188968863525798e+00, 1.1274396718809369e+01,  -1.8790661198015616e+01,
    2.6306925677221862e+01,  -3.0354145012179075e+01, 2.7463274058638209e+01,
    -1.6477964435182926e+01, 0.0000000000000000e+00,  1.6477964435182926e+01,
    -2.7463274058638209e+01, 3.0354145012179075e+01,  -2.6306925677221862e+01,
    1.8790661198015616e+01,  -1.1274396718809369e+01, 5.7188968863525798e+00,
    -2.4509558084368197e+00, 8.8234409103725508e-01,  -2.6394908706242670e-01,
    6.4520887948593189e-02,  -1.2569004145829843e-02, 1.8781270562734248e-03,
    -2.0225983682944576e-04, 1.3981094711712842e-05,  -4.6603649039042801e-07,
    2.5941169713712034e-06,  -7.2797407509104397e-05, 9.8025195055689352e-04,
    -8.4235842029744285e-03, 5.1814243856925575e-02,  -2.4246903280016285e-01,
    8.9516394110205377e-01,  -2.6664070668800939e+00, 6.4913561768485888e+00,
    -1.2978128062611834e+01, 2.1214265426542749e+01,  -2.7814727731222547e+01,
    2.7751464514244784e+01,  -1.7656926071870146e+01, -3.0820541604550561e-01,
    1.7845093589034775e+01,  -2.5149561949313256e+01, 1.7845093589034775e+01,
    -3.0820541604550561e-01, -1.7656926071870146e+01, 2.7751464514244784e+01,
    -2.7814727731222547e+01, 2.1214265426542749e+01,  -1.2978128062611834e+01,
    6.4913561768485888e+00,  -2.6664070668800939e+00, 8.9516394110205377e-01,
    -2.4246903280016285e-01, 5.1814243856925575e-02,  -8.4235842029744285e-03,
    9.8025195055689352e-04,  -7.2797407509104397e-05, 2.5941169713712034e-06,
    1.1691232321542031e-05,  -3.0616414642038194e-04, 3.8259557772246299e-03,
    -3.0311711899658006e-02, 1.7054000587433360e-01,  -7.2282432161771792e-01,
    2.3869667522003519e+00,  -6.2552031404138191e+00, 1.3090555918107396e+01,
    -2.1722554438601829e+01, 2.7751317407199895e+01,  -2.4995193418128675e+01,
    1.0644340923309549e+01,  9.5696718877831053e+00,  -2.2780059283676206e+01,
    1.8751878068879805e+01,  0.0000000000000000e+00,  -1.8751878068879805e+01,
    2.2780059283676206e+01,  -9.5696718877831053e+00, -1.0644340923309549e+01,
    2.4995193418128675e+01,  -2.7751317407199895e+01, 2.1722554438601829e+01,
    -1.3090555918107396e+01, 6.2552031404138191e+00,  -2.3869667522003519e+00,
    7.2282432161771792e-01,  -1.7054000587433360e-01, 3.0311711899658006e-02,
    -3.8259557772246299e-03, 3.0616414642038194e-04,  -1.1691232321542031e-05,
    4.5241558967149093e-05,  -1.1027629998242592e-03, 1.2747356515917244e-02,
    -9.2720933272403092e-02, 4.7449675410899761e-01,  -1.8076040501785993e+00,
    5.2805464198507774e+00,  -1.1969581694292478e+01, 2.0929303400664619e+01,
    -2.7287443554651361e+01, 2.3767281767212349e+01,  -7.4577015838115095e+00,
    -1.3413226523862413e+01, 2.2705033677997015e+01,  -1.1301506525198945e+01,
    -1.0816011694045546e+01, 2.1954889408808871e+01,  -1.0816011694045546e+01,
    -1.1301506525198945e+01, 2.2705033677997015e+01,  -1.3413226523862413e+01,
    -7.4577015838115095e+00, 2.3767281767212349e+01,  -2.7287443554651361e+01,
    2.0929303400664619e+01,  -1.1969581694292478e+01, 5.2805464198507774e+00,
    -1.8076040501785993e+00, 4.7449675410899761e-01,  -9.2720933272403092e-02,
    1.2747356515917244e-02,  -1.1027629998242592e-03, 4.5241558967149093e-05,
    1.5522490188229926e-04,  -3.5119634050870211e-03, 3.7418589956570554e-02,
    -2.4870096225170243e-01, 1.1501539322615948e+00,  -3.9016261036399991e+00,
    9.9402265465181259e+00,  -1.9029175320746088e+01, 2.6538951692463691e+01,
    -2.4154669033742763e+01, 7.6327558038613033e+00,  1.3999467890348040e+01,
    -2.1884870999228660e+01, 6.9693946747849482e+00,  1.5365058964657567e+01,
    -1.9385822058212817e+01, 0.0000000000000000e+00,  1.9385822058212817e+01,
    -1.5365058964657567e+01, -6.9693946747849482e+00, 2.1884870999228660e+01,
    -1.3999467890348040e+01, -7.6327558038613033e+00, 2.4154669033742763e+01,
    -2.6538951692463691e+01, 1.9029175320746088e+01,  -9.9402265465181259e+00,
    3.9016261036399991e+00,  -1.1501539322615948e+00, 2.4870096225170243e-01,
    -3.7418589956570554e-02, 3.5119634050870211e-03,  -1.5522490188229926e-04,
    4.8185679100586412e-04,  -1.0088876561685280e-02, 9.8687379551491344e-02,
    -5.9613942673884168e-01, 2.4724421681278552e+00,  -7.3825502391665569e+00,
    1.6091840692152612e+01,  -2.5081900159591495e+01, 2.5487465225487519e+01,
    -1.0523310289681163e+01, -1.2223152411073594e+01, 2.1618538930469395e+01,
    -6.1057667224516852e+00, -1.6491446666201796e+01, 1.6879676545447197e+01,
    5.8692537283905315e+00,  -2.0208063469901575e+01, 5.8692537283905315e+00,
    1.6879676545447197e+01,  -1.6491446666201796e+01, -6.1057667224516852e+00,
    2.1618538930469395e+01,  -1.2223152411073594e+01, -1.0523310289681163e+01,
    2.5487465225487519e+01,  -2.5081900159591495e+01, 1.6091840692152612e+01,
    -7.3825502391665569e+00, 2.4724421681278552e+00,  -5.9613942673884168e-01,
    9.8687379551491344e-02,  -1.0088876561685280e-02, 4.8185679100586412e-04,
    1.3722780701947185e-03,  -2.6502120230635501e-02, 2.3695590241914691e-01,
    -1.2928989772270643e+00, 4.7648371954404585e+00,  -1.2338354097465723e+01,
    2.2386851966397579e+01,  -2.6656941417187973e+01, 1.5348266209610093e+01,
    8.0550647703186744e+00,  -2.1558516149655606e+01, 8.2014503201049713e+00,
    1.5625931183651909e+01,  -1.6242861476963139e+01, -7.7016063331055165e+00,
    1.9254015832763791e+01,  0.0000000000000000e+00,  -1.9254015832763791e+01,
    7.7016063331055165e+00,  1.6242861476963139e+01,  -1.5625931183651909e+01,
    -8.2014503201049713e+00, 2.1558516149655606e+01,  -8.0550647703186744e+00,
    -1.5348266209610093e+01, 2.6656941417187973e+01,  -2.2386851966397579e+01,
    1.2338354097465723e+01,  -4.7648371954404585e+00, 1.2928989772270643e+00,
    -2.3695590241914691e-01, 2.6502120230635501e-02,  -1.3722780701947185e-03,
    3.6217967807219459e-03,  -6.4286892857814526e-02, 5.2279468208179114e-01,
    -2.5580692246161996e+00, 8.2829978716501831e+00,  -1.8232965984229885e+01,
    2.6381507168330124e+01,  -2.0878710690369619e+01, -1.2493768706775299e+00,
    2.0396423529067544e+01,  -1.2483407910478503e+01, -1.2782236360100359e+01,
    1.7364689557565221e+01,  6.5116572372290422e+00,  -1.8784518498202655e+01,
    -1.9625616341405758e+00, 1.9064884445937022e+01,  -1.9625616341405758e+00,
    -1.8784518498202655e+01, 6.5116572372290422e+00,  1.7364689557565221e+01,
    -1.2782236360100359e+01, -1.2483407910478503e+01, 2.0396423529067544e+01,
    -1.2493768706775299e+00, -2.0878710690369619e+01, 2.6381507168330124e+01,
    -1.8232965984229885e+01, 8.2829978716501831e+00,  -2.5580692246161996e+00,
    5.2279468208179114e-01,  -6.4286892857814526e-02, 3.6217967807219459e-03,
    8.9267041267638743e-03,  -1.4505894205991293e-01, 1.0669571117965111e+00,
    -4.6421021145654429e+00, 1.3018458548085839e+01,  -2.3678427749598690e+01,
    2.5290263831281820e+01,  -7.8267559421518200e+00, -1.6361397535637412e+01,
    1.7548616913319684e+01,  7.2778629342860608e+00,  -1.9084342664489306e+01,
    -2.5966307760192344e+00, 1.8896613778815670e+01,  6.9102271998418607e-01,
    -1.8657613439573019e+01, 0.0000000000000000e+00,  1.8657613439573019e+01,
    -6.9102271998418607e-01, -1.8896613778815670e+01, 2.5966307760192344e+00,
    1.9084342664489306e+01,  -7.2778629342860608e+00, -1.7548616913319684e+01,
    1.6361397535637412e+01,  7.8267559421518200e+00,  -2.5290263831281820e+01,
    2.3678427749598690e+01,  -1.3018458548085839e+01, 4.6421021145654429e+00,
    -1.0669571117965111e+00, 1.4505894205991293e-01,  -8.9267041267638743e-03,
    2.0670095438142064e-02,  -3.0617578867747930e-01, 2.0239190625984023e+00,
    -7.7498272140203399e+00, 1.8472374556657293e+01,  -2.6565097730985453e+01,
    1.7488564809808715e+01,  8.0791117210474273e+00,  -2.0923643243975334e+01,
    1.2375256840198827e+00,  1.9327197389815126e+01,  -3.8915521428854851e+00,
    -1.8374933105637439e+01, 3.3477751591281439e+00,  1.8212173348301871e+01,
    -1.2746036736260380e+00, -1.8246957854014862e+01, -1.2746036736260380e+00,
    1.8212173348301871e+01,  3.3477751591281439e+00,  -1.8374933105637439e+01,
    -3.8915521428854851e+00, 1.9327197389815126e+01,  1.2375256840198827e+00,
    -2.0923643243975334e+01, 8.0791117210474273e+00,  1.7488564809808715e+01,
    -2.6565097730985453e+01, 1.8472374556657293e+01,  -7.7498272140203399e+00,
    2.0239190625984023e+00,  -3.0617578867747930e-01, 2.0670095438142064e-02,
    4.5180960042357270e-02,  -6.0711915056917587e-01, 3.5802267207758272e+00,
    -1.1913654401169083e+01, 2.3514440705426384e+01,  -2.4728506248444443e+01,
    4.1020944505643007e+00,  1.9351488280667098e+01,  -1.1473501262280354e+01,
    -1.5474613879946869e+01, 1.1851609107529159e+01,  1.5197013568029222e+01,
    -9.4515276581155589e+00, -1.6547014623265930e+01, 5.2570987038384445e+00,
    1.7742708125454751e+01,  0.0000000000000000e+00,  -1.7742708125454751e+01,
    -5.2570987038384445e+00, 1.6547014623265930e+01,  9.4515276581155589e+00,
    -1.5197013568029222e+01, -1.1851609107529159e+01, 1.5474613879946869e+01,
    1.1473501262280354e+01,  -1.9351488280667098e+01, -4.1020944505643007e+00,
    2.4728506248444443e+01,  -2.3514440705426384e+01, 1.1913654401169083e+01,
    -3.5802267207758272e+00, 6.0711915056917587e-01,  -4.5180960042357270e-02,
    9.3587837517216846e-02,  -1.1347525298962542e+00, 5.9182986120262564e+00,
    -1.6841659679660768e+01, 2.6478174968664145e+01,  -1.7062772657287653e+01,
    -1.0453934086747147e+01, 1.9553121764139465e+01,  5.7979868837534045e+00,
    -1.8212112501680124e+01, -7.5038264866713851e+00, 1.5696640074708744e+01,
    1.1815386409932231e+01,  -1.1032837496587840e+01, -1.5945418200077786e+01,
    4.0187639366049703e+00,  1.7630706302525031e+01,  4.0187639366049703e+00,
    -1.5945418200077786e+01, -1.1032837496587840e+01, 1.1815386409932231e+01,
    1.5696640074708744e+01,  -7.5038264866713851e+00, -1.8212112501680124e+01,
    5.7979868837534045e+00,  1.9553121764139465e+01,  -1.0453934086747147e+01,
    -1.7062772657287653e+01, 2.6478174968664145e+01,  -1.6841659679660768e+01,
    5.9182986120262564e+00,  -1.1347525298962542e+00, 9.3587837517216846e-02,
    1.8430030171174749e-01,  -2.0042657811152540e+00, 9.1503613507126893e+00,
    -2.1800198994814970e+01, 2.5600880202876404e+01,  -4.6120432982724102e+00,
    -1.9903604019366966e+01, 8.0422019726227543e+00,  1.8192701677350168e+01,
    -4.4143663617590043e+00, -1.8352051052972740e+01, -2.7813406349955421e+00,
    1.6365569533452028e+01,  1.1004434686286709e+01,  -9.9563932875927357e+00,
    -1.6593988812654562e+01, 0.0000000000000000e+00,  1.6593988812654562e+01,
    9.9563932875927357e+00,  -1.1004434686286709e+01, -1.6365569533452028e+01,
    2.7813406349955421e+00,  1.8352051052972740e+01,  4.4143663617590043e+00,
    -1.8192701677350168e+01, -8.0422019726227543e+00, 1.9903604019366966e+01,
    4.6120432982724102e+00,  -2.5600880202876404e+01, 2.1800198994814970e+01,
    -9.1503613507126893e+00, 2.0042657811152540e+00,  -1.8430030171174749e-01,
    3.4596507459410830e-01,  -3.3515366601304244e+00, 1.3227584021376350e+01,
    -2.5616063234815780e+01, 1.9769845155267888e+01,  9.0804046563256602e+00,
    -1.932712659253413e+01,  -8.2949439973067314e+00, 1.5928696399744720e+01,
    1.3137557097418776e+01,  -9.0186629748404634e+00, -1.7330844293948672e+01,
    -1.7927456062254210e+00, 1.5625407726987476e+01,  1.2560860448813177e+01,
    -6.3687093641070831e+00, -1.7151375713238878e+01, -6.3687093641070831e+00,
    1.2560860448813177e+01,  1.5625407726987476e+01,  -1.7927456062254210e+00,
    -1.7330844293948672e+01, -9.0186629748404634e+00, 1.3137557097418776e+01,
    1.5928696399744720e+01,  -8.2949439973067314e+00, -1.9327126592534139e+01,
    9.0804046563256602e+00,  1.9769845155267888e+01,  -2.5616063234815780e+01,
    1.3227584021376350e+01,  -3.3515366601304244e+00, 3.4596507459410830e-01,
};

template <typename T> struct CompareSecondElementOfPair
{
    bool operator()( const T& inLeft, const T& inRight )
    {
        return inLeft.second < inRight.second;
    }
};

template <std::size_t MaxNInterval_> class Errors
{
private:
    std::array<double, MaxNInterval_> mData;
    std::priority_queue<
        std::pair<std::size_t, double>,
        std::vector<std::pair<std::size_t, double>>,
        CompareSecondElementOfPair<std::pair<std::size_t, double>>>
        mQue;

public:
    Errors() {}
    void set( std::size_t inInd, double inVal )
    {
        mQue.push( { inInd, inVal } );
        mData[inInd] = inVal;
    }
    std::pair<std::size_t, double> maxError()
    {
        // auto lItr = std::max_element( mData.begin(), mData.end() );
        // return { lItr - mData.begin(), *lItr };
        std::pair<std::size_t, double> lResult = mQue.top();
        while ( !mQue.empty() )
        {
            lResult = mQue.top();
            if ( mData[lResult.first] == lResult.second ) { break; }
            mQue.pop();
        }
        return lResult;
    }
};

template <typename Func_, typename TypeY_, std::size_t MaxNInterval_,
          std::size_t MaxNFuncEval_>
class DoublyAdaptiveNewtonCotes
{
private:
    Func_ mFunc;
    std::size_t mNEvalFunc;
    std::size_t mSizeInterval;

    std::vector<short> mFlags;
    std::vector<double> mLeftPos;
    std::vector<double> mRightPos;
    std::vector<TypeY_> mIntegral;
    std::vector<double> mAbsIntegral;
    Errors<MaxNInterval_> mErrors;
    std::vector<std::vector<TypeY_>> mVals;

    double mTotalError;
    TypeY_ mTotalIntegral;
    double mTotalAbsIntegral;

    std::vector<TypeY_> applyFunc( const std::vector<double>& inXs )
    {
        std::vector<TypeY_> lResult( inXs.size() );
        for ( std::size_t iX = 0; iX < inXs.size(); ++iX )
        {
            lResult[iX] = mFunc( inXs[iX] );
        }
        mNEvalFunc += inXs.size();
        return lResult;
    }
    double getNoise() const { return gFactorEps * mTotalAbsIntegral; }
    void subtractTotal( std::size_t inIndInterval, double inErr )
    {
        assert( inIndInterval < MaxNInterval_ );
        mTotalError -= inErr;
        mTotalIntegral -= mIntegral[inIndInterval];
        mTotalAbsIntegral -= mAbsIntegral[inIndInterval];
    }
    void updateElement( std::size_t inIndInterval, double inErr, double inInt,
                        double inAbsInt )
    {
        assert( inIndInterval < MaxNInterval_ );
        mIntegral[inIndInterval]    = inInt;
        mAbsIntegral[inIndInterval] = inAbsInt;
        addTotal( inIndInterval, inErr );
    }
    void addTotal( std::size_t inIndInterval, double inErr )
    {
        assert( inIndInterval < MaxNInterval_ );
        mErrors.set( inIndInterval, inErr );
        mTotalError += inErr;
        mTotalIntegral += mIntegral[inIndInterval];
        mTotalAbsIntegral += mAbsIntegral[inIndInterval];
    }
    template <std::size_t NEach_>
    std::pair<std::vector<TypeY_>, std::vector<TypeY_>>
    previousValsForBisectCotes( std::size_t inIndInterval )
    {
        assert( inIndInterval < MaxNInterval_ );
        std::vector<double> lLeftVals( NEach_ ), lRightVals( NEach_ );
        constexpr std::size_t lIntervalRef = ( ( 16 + NEach_ ) / NEach_ );
        for ( std::size_t iVal = 0, iRef = 0; iVal < NEach_;
              ++iVal, iRef += lIntervalRef )
        {
            assert( iRef + 16 < 33 );
            lLeftVals[iVal]  = mVals[inIndInterval][iRef];
            lRightVals[iVal] = mVals[inIndInterval][iRef + 16];
        }
        return { lLeftVals, lRightVals };
    }
    template <std::size_t NPoints_>
    std::pair<std::tuple<TypeY_, double, double, double>,
              std::tuple<TypeY_, double, double, double>>
    bisectCotes( std::size_t inIndInterval,
                 const std::vector<TypeY_>& inLeftVals,
                 const std::vector<TypeY_>& inRightVals )
    {
        assert( inIndInterval < MaxNInterval_ );
        assert( inLeftVals.size() == NPoints_ &&
                inRightVals.size() == NPoints_ );
        const double lTmpMid =
            0.5 * ( mLeftPos[inIndInterval] + mRightPos[inIndInterval] );
        return { evalNewtonCotesIntegralAndErr<NPoints_, TypeY_>(
                     mLeftPos[inIndInterval], lTmpMid, inLeftVals ),
                 evalNewtonCotesIntegralAndErr<NPoints_, TypeY_>(
                     lTmpMid, mRightPos[inIndInterval], inRightVals ) };
    }
    std::tuple<TypeY_, double, double, double> intervalCotes9(
        std::size_t inIndInterval )
    {
        assert( inIndInterval < MaxNInterval_ );
        std::vector<double> lTmpX( 4 );
        std::vector<TypeY_> lCotesVals( 9 );
        const double lDif =
            0.25 * ( mRightPos[inIndInterval] - mLeftPos[inIndInterval] );
        lTmpX[0] = mLeftPos[inIndInterval] + 0.5 * lDif;
        for ( std::size_t iX = 0; iX < 3; ++iX )
        {
            lTmpX[iX + 1] = lTmpX[iX] + lDif;
        }
        auto lTmpY = applyFunc( lTmpX );
        for ( std::size_t iY = 0, iRef = 4; iY < 4; ++iY, iRef += 8 )
        {
            mVals[inIndInterval][iRef] = lTmpY[iY];
        }
        for ( std::size_t iVal = 0, iRef = 0; iVal < 9; ++iVal, iRef += 4 )
        {
            lCotesVals[iVal] = mVals[inIndInterval][iRef];
        }
        return evalNewtonCotesIntegralAndErr<9, TypeY_>(
            mLeftPos[inIndInterval], mRightPos[inIndInterval], lCotesVals );
    }
    std::tuple<TypeY_, double, double, double> intervalCotes17(
        std::size_t inIndInterval )
    {
        assert( inIndInterval < MaxNInterval_ );
        std::vector<double> lTmpX( 8 );
        std::vector<TypeY_> lCotesVals( 17 );
        const double lDif =
            0.125 * ( mRightPos[inIndInterval] - mLeftPos[inIndInterval] );
        lTmpX[0] = mLeftPos[inIndInterval] + 0.5 * lDif;
        for ( std::size_t iX = 0; iX < 7; ++iX )
        {
            lTmpX[iX + 1] = lTmpX[iX] + lDif;
        }
        auto lTmpY = applyFunc( lTmpX );
        for ( std::size_t iY = 0, iRef = 2; iY < 8; ++iY, iRef += 4 )
        {
            mVals[inIndInterval][iRef] = lTmpY[iY];
        }
        for ( std::size_t iVal = 0, iRef = 0; iVal < 17; ++iVal, iRef += 2 )
        {
            lCotesVals[iVal] = mVals[inIndInterval][iRef];
        }
        return evalNewtonCotesIntegralAndErr<17, TypeY_>(
            mLeftPos[inIndInterval], mRightPos[inIndInterval], lCotesVals );
    }
    std::tuple<TypeY_, double, double, double> intervalCotes33(
        std::size_t inIndInterval )
    {
        assert( inIndInterval < MaxNInterval_ );
        std::vector<double> lTmpX( 16 );
        std::vector<TypeY_> lCotesVals( 33 );
        const double lDif =
            0.0625 * ( mRightPos[inIndInterval] - mLeftPos[inIndInterval] );
        lTmpX[0] = mLeftPos[inIndInterval] + 0.5 * lDif;
        for ( std::size_t iX = 0; iX < 15; ++iX )
        {
            lTmpX[iX + 1] = lTmpX[iX] + lDif;
        }
        auto lTmpY = applyFunc( lTmpX );
        for ( std::size_t iY = 0, iRef = 1; iY < 16; ++iY, iRef += 2 )
        {
            mVals[inIndInterval][iRef] = lTmpY[iY];
        }
        return evalNewtonCotesIntegralAndErr<33, TypeY_>(
            mLeftPos[inIndInterval], mRightPos[inIndInterval],
            mVals[inIndInterval] );
    }
    template <std::size_t NEach_, short Flag_>
    void addInterval( std::size_t inIndInterval, double inLErr, double inLInt,
                      double inLAbsInt, const std::vector<TypeY_>& inLeftVals,
                      double inRErr, double inRInt, double inRAbsInt,
                      const std::vector<TypeY_>& inRightVals )
    {
        assert( inIndInterval < MaxNInterval_ );
        updateElement( inIndInterval, inLErr, inLInt, inLAbsInt );
        updateElement( mSizeInterval, inRErr, inRInt, inRAbsInt );

        mVals[inIndInterval].assign( 33, 0.0 );
        for ( std::size_t iVal = 0, iRef = 0; iVal < NEach_;
              ++iVal, iRef += ( 32 / ( NEach_ - 1 ) ) )
        {
            assert( iRef < 33 );
            mVals[inIndInterval][iRef] = inLeftVals[iVal];
            mVals[mSizeInterval][iRef] = inRightVals[iVal];
        }

        mRightPos[mSizeInterval] = mRightPos[inIndInterval];
        mRightPos[inIndInterval] =
            0.5 * ( mLeftPos[inIndInterval] + mRightPos[inIndInterval] );
        mLeftPos[mSizeInterval] = mRightPos[inIndInterval];

        mFlags[inIndInterval] = Flag_;
        mFlags[mSizeInterval] = Flag_;

        ++mSizeInterval;
    }
    void step()
    {
        auto [lIndInterval, lErr] = mErrors.maxError();
        subtractTotal( lIndInterval, lErr );
        const short lTmpFlag = mFlags[lIndInterval];
        if ( lTmpFlag == 1 || lTmpFlag == 2 )
        {
            const auto [lLeftVals, lRightVals] =
                previousValsForBisectCotes<5>( lIndInterval );
            const auto [lLeftTuple, lRightTuple] =
                bisectCotes<5>( lIndInterval, lLeftVals, lRightVals );
            auto& [lLInt, lLAbsInt, lLErr, lLRMax] = lLeftTuple;
            auto& [lRInt, lRAbsInt, lRErr, lRRMax] = lRightTuple;
            if ( lTmpFlag == 2 && ( lLErr + lRErr >= lErr ) )
            {
                double lRMax;
                std::tie( mIntegral[lIndInterval], mAbsIntegral[lIndInterval],
                          lErr, lRMax ) = intervalCotes17( lIndInterval );
                addTotal( lIndInterval, lErr );
                mFlags[lIndInterval] = 4;
            }
            else
            {
                addInterval<5, 0>( lIndInterval, lLErr, lLInt, lLAbsInt,
                                   lLeftVals, lRErr, lRInt, lRAbsInt,
                                   lRightVals );
            }
        }
        else if ( lTmpFlag == 0 )
        {
            double lRMax;
            std::tie( mIntegral[lIndInterval], mAbsIntegral[lIndInterval], lErr,
                      lRMax ) = intervalCotes9( lIndInterval );
            addTotal( lIndInterval, lErr );
            mFlags[lIndInterval] = 1 + ( lRMax < 0.25 );
        }
        else if ( lTmpFlag == 3 )
        {
            double lRMax;
            std::tie( mIntegral[lIndInterval], mAbsIntegral[lIndInterval], lErr,
                      lRMax ) = intervalCotes17( lIndInterval );
            addTotal( lIndInterval, lErr );
            mFlags[lIndInterval] = 4 + ( lRMax < 0.125 );
        }
        else if ( lTmpFlag == 4 || lTmpFlag == 5 )
        {
            auto [lLeftVals, lRightVals] =
                previousValsForBisectCotes<9>( lIndInterval );
            auto [lLeftTuple, lRightTuple] =
                bisectCotes<9>( lIndInterval, lLeftVals, lRightVals );
            auto& [lLInt, lLAbsInt, lLErr, lLRMax] = lLeftTuple;
            auto& [lRInt, lRAbsInt, lRErr, lRRMax] = lRightTuple;
            if ( lTmpFlag == 5 && ( lLErr + lRErr >= lErr ) )
            {
                std::tie( mIntegral[lIndInterval], mAbsIntegral[lIndInterval],
                          lErr, std::ignore ) = intervalCotes33( lIndInterval );
                addTotal( lIndInterval, lErr );
                mFlags[lIndInterval] = 7;
            }
            else
            {
                addInterval<9, 3>( lIndInterval, lLErr, lLInt, lLAbsInt,
                                   lLeftVals, lRErr, lRInt, lRAbsInt,
                                   lRightVals );
            }
        }
        else if ( lTmpFlag == 6 )
        {
            std::tie( mIntegral[lIndInterval], mAbsIntegral[lIndInterval], lErr,
                      std::ignore ) = intervalCotes33( lIndInterval );
            addTotal( lIndInterval, lErr );
            mFlags[lIndInterval] = 7;
        }
        else if ( lTmpFlag == 7 )
        {
            auto [lLeftVals, lRightVals] =
                previousValsForBisectCotes<17>( lIndInterval );
            auto [lLeftTuple, lRightTuple] =
                bisectCotes<17>( lIndInterval, lLeftVals, lRightVals );
            auto& [lLInt, lLAbsInt, lLErr, lLRMax] = lLeftTuple;
            auto& [lRInt, lRAbsInt, lRErr, lRRMax] = lRightTuple;
            addInterval<17, 6>( lIndInterval, lLErr, lLInt, lLAbsInt, lLeftVals,
                                lRErr, lRInt, lRAbsInt, lRightVals );
        }
    }

public:
    ~DoublyAdaptiveNewtonCotes() {}
    DoublyAdaptiveNewtonCotes( Func_ inFunc, double inMin, double inMax ) :
        mFunc( inFunc ),
        mNEvalFunc( 0 ),
        mSizeInterval( 1 ),
        mFlags( MaxNInterval_, 0 ),
        mLeftPos( MaxNInterval_ ),
        mRightPos( MaxNInterval_ ),
        mIntegral( MaxNInterval_ ),
        mAbsIntegral( MaxNInterval_ ),
        mErrors(),
        mVals( MaxNInterval_, std::vector<double>( 33, 0.0 ) )
    {
        std::vector<double> lTmpX( 9 );
        const double lTmpDif = 0.125 * ( inMax - inMin );
        lTmpX[0]             = inMin;
        for ( std::size_t j = 0; j < 8; ++j )
        {
            lTmpX[j + 1] = lTmpX[j] + lTmpDif;
        }
        const auto lTmpY = applyFunc( lTmpX );

        double lTmpErr, lTmpRMax;
        std::tie( mIntegral[0], mAbsIntegral[0], lTmpErr, lTmpRMax ) =
            evalNewtonCotesIntegralAndErr<9, TypeY_>( inMin, inMax, lTmpY );

        mErrors.set( 0, lTmpErr );
        mTotalError       = lTmpErr;
        mTotalIntegral    = mIntegral[0];
        mTotalAbsIntegral = mAbsIntegral[0];
        for ( std::size_t iY = 0, iVal = 0; iY < 9; ++iY, iVal += 4 )
        {
            mVals[0][iVal] = lTmpY[iY];
        }
        mLeftPos[0]  = inMin;
        mRightPos[0] = inMax;
        mFlags[0]    = 1 + ( lTmpRMax < 0.25 );
    }
    double operator()( double inTolAbs, double inTolRel )
    {
        while (
            ( mTotalError > std::max( { std::abs( mTotalIntegral ) * inTolRel,
                                        inTolAbs, getNoise() } ) ) &&
            ( mNEvalFunc < MaxNFuncEval_ - 3 ) &&
            ( mSizeInterval < MaxNInterval_ - 3 ) )
        {
            step();
        }
        return mTotalIntegral;
    }
};

template <typename TypeY_>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inMin, double inMax,
                                  double inTolAbs, double inTolRel )
{
    DoublyAdaptiveNewtonCotes<decltype( inFunc ), TypeY_> lObj( inFunc, inMin,
                                                                inMax );
    return lObj( inTolAbs, inTolRel );
}

template <std::size_t N_, std::array<double, N_> const& Weight_,
          typename TypeY_>
double weightedSum( const std::vector<TypeY_>& inVals )
{
    assert( values.size() == N_ );
    double lSum = 0.0;
    for ( std::size_t i = 0; i < N_; ++i ) { lSum += Weight_[i] * inVals[i]; }
    return lSum;
}

template <std::size_t N_, std::array<double, N_> const& Weight_,
          typename TypeY_>
double weightedAbsSum( const std::vector<TypeY_>& inVals )
{
    assert( values.size() == N_ );
    double lSum = 0.0;
    for ( std::size_t i = 0; i < N_; ++i )
    {
        lSum += std::abs( Weight_[i] * inVals[i] );
    }
    return lSum;
}

template <std::size_t N1_, std::size_t N2_,
          std::array<double, N1_ * N2_> const& Weight_, typename TypeY_>
std::vector<double> weightedSum2D( const std::vector<TypeY_>& inVals )
{
    assert( values.size() == N2_ );
    std::vector<double> lResult( N1_ );
    for ( std::size_t i = 0, iN2 = 0; i < N1_; ++i, iN2 += N2_ )
    {
        for ( std::size_t j = 0; j < N2_; ++j )
        {
            lResult[i] += Weight_[iN2 + j] * inVals[j];
        }
    }
    return lResult;
}

template <std::size_t N_, std::array<double, N_> const& Weight_,
          typename TypeY_>
std::tuple<double, double> evalIntAndAbsInt( double lHalfDif,
                                             const std::vector<TypeY_>& inVals )
{
    return {
        lHalfDif * weightedSum<N_, Weight_>( inVals ),
        lHalfDif * weightedAbsSum<N_, Weight_>( inVals ),
    };
}

template <std::size_t N1_, std::size_t N2_,
          std::array<double, N1_ * N2_> const& Weight_, double const& Coeff_,
          typename TypeY_>
std::pair<double, double> evalErrAndRMax( double inHalfDif, double inAbsInt,
                                          const std::vector<TypeY_>& inVals )
{
    std::vector<double> lTmpErrVec = weightedSum2D<N1_, N2_, Weight_>( inVals );

    std::size_t lNErrVec;
    if constexpr ( N1_ < 5 ) { lNErrVec = N1_; }
    else if constexpr ( N1_ < 9 ) { lNErrVec = N1_ / 2; }
    else { lNErrVec = 5; }

    std::vector<double> lErrVec( lNErrVec );
    if constexpr ( N1_ < 5 )
    {
        for ( std::size_t iRes = 0; iRes < N1_; ++iRes )
        {
            lErrVec[iRes] = inHalfDif * std::abs( lTmpErrVec[iRes] );
        }
    }
    else if constexpr ( N1_ < 9 )
    {
        for ( std::size_t iRes = 0, iE = 0; iRes < lNErrVec; ++iRes, iE += 2 )
        {
            lErrVec[iRes] =
                inHalfDif *
                std::sqrt( lTmpErrVec[iE] * lTmpErrVec[iE] +
                           lTmpErrVec[iE + 1] * lTmpErrVec[iE + 1] );
        }
    }
    else
    {
        for ( std::size_t iRes = 0, iE = 0; iRes < 5; ++iRes, iE += 3 )
        {
            lErrVec[iRes] =
                inHalfDif *
                std::sqrt( lTmpErrVec[iE] * lTmpErrVec[iE] +
                           lTmpErrVec[iE + 1] * lTmpErrVec[iE + 1] +
                           lTmpErrVec[iE + 2] * lTmpErrVec[iE + 2] );
        }
    }

    double lRMax = 0.0;
    for ( std::size_t i = 0; i < lNErrVec - 1; ++i )
    {
        if ( lErrVec[i + 1] == 0.0 )
        {
            lRMax = 2.0;
            break;
        }
        lRMax = std::max( { lRMax, lErrVec[i] / lErrVec[i + 1] } );
    }
    const double lNoise = gFactorEps * inAbsInt;
    if ( lErrVec[0] < lNoise && lErrVec[1] < lNoise ) { return { 0.0, lRMax }; }
    if ( inHalfDif <= gFactorEps ) { return { 0.0, lRMax }; }
    if ( lRMax > 1.0 )
    {
        return {
            Coeff_ * ( *std::max_element( lErrVec.begin(), lErrVec.end() ) ),
            lRMax };
    }

    if constexpr ( N1_ < 5 )
    {
        if ( lRMax > 0.5 ) { return { Coeff_ * lRMax * lErrVec[1], lRMax }; }
        return { Coeff_ * std::pow( 2 * lRMax, 3 ) * lRMax * lErrVec[1],
                 lRMax };
    }
    else if constexpr ( N1_ < 9 )
    {
        if ( lRMax > 0.25 ) { return { Coeff_ * lRMax * lErrVec[0], lRMax }; }
        return { 4.0 * Coeff_ * lRMax * lRMax * lErrVec[0], lRMax };
    }
    else if constexpr ( N2_ < 33 )
    {
        if ( lRMax > 0.125 )
        {
            return { 4.0 * Coeff_ * lRMax * lRMax * lErrVec[0], lRMax };
        }
        return { Coeff_ * pow( 8.0 * lRMax, 2.0 / 3.0 ) * lRMax * lErrVec[0],
                 lRMax };
    }
    return { Coeff_ * inHalfDif *
                 std::sqrt( lTmpErrVec[0] * lTmpErrVec[0] +
                            lTmpErrVec[1] * lTmpErrVec[1] +
                            lTmpErrVec[2] * lTmpErrVec[2] +
                            lTmpErrVec[3] * lTmpErrVec[3] ),
             lRMax };
}

template <std::size_t NPoints_, typename TypeY_>
std::tuple<TypeY_, double, double, double> evalNewtonCotesIntegralAndErr(
    double inLeft, double inRight, const std::vector<TypeY_>& inVals )
{
    assert( inVals.size() == NPoints_ );
    const double lHalfDif = 0.5 * ( inRight - inLeft );
    if constexpr ( NPoints_ == 5 )
    {
        const auto [lIntegral, lAbsIntegral] =
            evalIntAndAbsInt<NPoints_, gWeightForCotes5>( lHalfDif, inVals );
        const auto [lErr, lRMax] =
            evalErrAndRMax<4, NPoints_, gNullWeightForCotes5, gCoeffD>(
                lHalfDif, lAbsIntegral, inVals );
        return { lIntegral, lAbsIntegral, lErr, lRMax };
    }
    else if constexpr ( NPoints_ == 9 )
    {
        const auto [lIntegral, lAbsIntegral] =
            evalIntAndAbsInt<NPoints_, gWeightForCotes9>( lHalfDif, inVals );
        const auto [lErr, lRMax] =
            evalErrAndRMax<8, NPoints_, gNullWeightForCotes9, gCoeffD>(
                lHalfDif, lAbsIntegral, inVals );
        return { lIntegral, lAbsIntegral, lErr, lRMax };
    }
    else if constexpr ( NPoints_ == 17 )
    {
        const auto [lIntegral, lAbsIntegral] =
            evalIntAndAbsInt<NPoints_, gWeightForCotes17>( lHalfDif, inVals );
        const auto [lErr, lRMax] =
            evalErrAndRMax<15, NPoints_, gNullWeightForCotes17, gCoeffD>(
                lHalfDif, lAbsIntegral, inVals );
        return { lIntegral, lAbsIntegral, lErr, lRMax };
    }
    else
    {
        const auto [lIntegral, lAbsIntegral] =
            evalIntAndAbsInt<NPoints_, gWeightForCotes33>( lHalfDif, inVals );
        const auto [lErr, lRMax] =
            evalErrAndRMax<15, NPoints_, gNullWeightForCotes33, gCoeffD>(
                lHalfDif, lAbsIntegral, inVals );
        return { lIntegral, lAbsIntegral, lErr, lRMax };
    }
}
}  // namespace FiniteInterval

namespace InfiniteInterval
{
template <typename TypeY_>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inTolAbs,
                                  double inTolRel )
{
    auto lTransfFunc = [inFunc]( double inX ) -> TypeY_
    {
        const double lInvX = 1.0 / inX;
        const double lArg  = ( 1.0 - inX ) * lInvX;
        return ( inFunc( lArg ) + inFunc( -lArg ) ) * lInvX * lInvX;
    };
    return FiniteInterval::doublyAdaptiveNewtonCotes(
        lTransfFunc, std::numeric_limits<double>::epsilon(), 1.0, inTolAbs,
        inTolRel );
}
}  // namespace InfiniteInterval

namespace UpperInfiniteInterval
{
template <typename TypeY_>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inMin, double inTolAbs,
                                  double inTolRel )
{
    auto lTransfFunc = [inFunc, inMin]( double inX ) -> TypeY_
    {
        const double lInvX = 1.0 / inX;
        return inFunc( inMin + ( 1.0 - inX ) * lInvX ) * lInvX * lInvX;
    };
    return FiniteInterval::doublyAdaptiveNewtonCotes(
        lTransfFunc, std::numeric_limits<double>::epsilon(), 1.0, inTolAbs,
        inTolRel );
}

}  // namespace UpperInfiniteInterval

}  // namespace Math::Integral
