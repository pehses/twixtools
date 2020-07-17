# remove fidnav:
twixzip.py -i meas_MID00359_FID50319_T2.dat -o /tmp/test_t2.h5 --testmode --remove_fidnav

# # remove_os:
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --remove_os

# zfp:
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-8
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-7
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-6
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-4

# # zfp prec:
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 6
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 7
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 8
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 9
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 10
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 12
# twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_prec 16


# # scc:
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --scc --ncc 8
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --scc --ncc 16
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --scc --ncc 24
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --scc --ncc 32

# # gcc:
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 4
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 6
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 8
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 16
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 24
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --gcc --ncc 32


# zfp & remove_os:
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-6 --remove_os

# zfp & remove_os & gcc
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-6 --remove_os --gcc --ncc 16

# zfp & remove_os & scc
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-6 --remove_os --scc --ncc 24



twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5 --remove_os --gcc --ncc 8
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5 --remove_os --gcc --ncc 16


twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5 --remove_os --scc --ncc 16
twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5 --remove_os --scc --ncc 24

twixzip.py -i meas_MID00359_FID50319_T2_rmfidnav.dat -o /tmp/test_t2.h5 --testmode --zfp --zfp_tol 1e-5 --remove_os

# results:
# original size: 7567681536
# rm_os: compressed size = 3819262352, comp. factor = 1.9814510862384453
# zfp_tol 1e-8: compressed size = 3820073460  comp. factor = 1.9810303689814384
# zfp_tol 1e-7: compressed size = 3116307567  comp. factor = 2.428412912813108
# zfp_tol 2e-7: compressed size = 2881703903  comp. factor = 2.626113504625392
# zfp_tol 5e-7: compressed size = 2412457486  comp. factor = 3.136918092823129
# zfp_tol 1e-6: compressed size = 2177810580  comp. factor = 3.4749034674999146
# zfp_tol 1e-5: compressed size = 1477185138  comp. factor = 5.123042021832202
# zfp_prec 6: compressed size = 1660357650  comp. factor = 4.557862299125733
# zfp_prec 7: compressed size = 1895936730  comp. factor = 3.991526413436803
# zfp_prec 8: compressed size = 2123538610  comp. factor = 3.5637127106438626
# zfp_prec 9: compressed size = 2358160138  comp. factor = 3.2091465774747143
# zfp_prec 10: compressed size = 2592787840  comp. factor = 2.918743068464869
# zfp_prec 12: compressed size = 3062016962  comp. factor = 2.4714695019380497
# zfp_prec 16: compressed size = 4000509277  comp. factor = 1.891679536780137
# scc 8: 
# scc 16: 
# scc 24: 
# scc 32: 
# gcc 8: 
# gcc 16: 
# gcc 24: 
# gcc 32: 
# zfp_tol 2e-7 + rm_os: compressed size = 961811141  comp. factor = 7.868157493093543
# zfp_tol 2e-7 + rm_os + scc 24: 
# zfp_tol 2e-7 + rm_os + gcc 16: 
