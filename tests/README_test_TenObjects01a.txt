INVOCATION FOR make_tracklets:

time make_tracklets_new -dets test_TenObjects01a.csv -outimgs imgs_test_TenObjects01a.txt -pairdets pairdets_test_TenObjects01a.csv -tracklets tracklets_test_TenObjects01a.csv -trk2det trk2det_test_TenObjects01a.csv -colformat colformat_LSST_02.txt -imrad 2.0 -maxtime 1.5 -maxGCR 0.5 -maxvel 1.5 -earth Earth1day2020s_02a.csv -obscode ObsCodesNew.txt

EXPECTED FINAL OUTPUT AND TIMING:

Writing paired detection file with 107 lines
Writing tracklet file with 41 lines
Writing trk2det file with 107 lines
0.033u 0.006s 0:00.08 37.5%	0+0k 0+72io 0pf+0w


INVOCATION FOR heliolinc:

time heliolinc_kd -imgs imgs_test_TenObjects01a.txt -pairdets pairdets_test_TenObjects01a.csv -tracklets tracklets_test_TenObjects01a.csv -trk2det trk2det_test_TenObjects01a.csv -mjd 61109.09 -obspos Earth1day2020s_02a.csv -heliodist heliohyp_rmb00a.txt -clustrad 200000.0 -outsum sum_test_TenObjects01a.csv -clust2det clust2det_test_TenObjects01a.csv

EXPECTED FINAL OUTPUT AND TIMING:

De-duplicating output set of 52 candidate linkages
Final de-duplicating set contains 25 linkages
Writing 25 lines to output cluster-summary file sum_test_TenObjects01a.csv
Writing 255 lines to output clust2det file clust2det_test_TenObjects01a.csv
0.008u 0.001s 0:00.02 0.0%	0+0k 0+24io 0pf+0w


SETUP AND INVOCATION FOR link_planarity:

printf "sum_test_TenObjects01a.csv clust2det_test_TenObjects01a.csv\n" > clusterlist_test_TenObjects01a

time link_planarity -imgs imgs_test_TenObjects01a.txt -pairdet pairdets_test_TenObjects01a.csv -lflist clusterlist_test_TenObjects01a -mjd 61109.09 -simptype 1 -max_astrom_rms 0.5 -oop 10000.0 -ptpow 3 -maxrms 400000.0 -outsum LPLsum_test_TenObjects01a.csv -clust2det LPLclust2det_test_TenObjects01a.csv

EXPECTED FINAL OUTPUT AND TIMING:

Writing 10 lines to output cluster-summary file LPLsum_test_TenObjects01a.csv
Writing 105 lines to output clust2det file LPLclust2det_test_TenObjects01a.csv
0.036u 0.000s 0:00.05 60.0%	0+0k 32+16io 0pf+0w
