import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
from matplotlib import rcParams
import numpy as np
import csv
import pathlib

isositei_list = ["'[23Na]+'","'[25Mg]+'","'[27Al]+'","'[28Si]+'","'[31P]+'","'[32S]+'","'[44Ca]+'","'[48Ti]+'","'[51V]+'","'[53Cr]+'","'[55Mn]+'","'[57Fe]+'","'[59Co]+'","'[60Ni]+'","'[195Pt]+'"]
# measure time
mea_time = 3  # [min]
# integration time
int_time = 0.002  # [sec]
# averaging N points
N = 300  # 移動平均の計算点数
# threshold level
n = 3  # 移動標準偏差からのしきいレベル
n2 = 3  # 粒子判別のしきいレベル
sum_thres = 0  # 全同位体合計値からの粒子判別
# end cycles
end = 2  # 移動平均・移動標準偏差の計算サイクル数
# number of bins
#bin = 20
date = "20211023"
sample = "MQ"
label1 = ["_1","_2","_3"]
datasheet2 = ["C:/python/{}/csv/{}{}.csv".format(date,sample, i) for i in label1]
csv_title1 = ["{}{}_subtract_ma.csv".format(sample, i) for i in label1]
csv_title5 = ["{}{}_NP_events_large.csv".format(sample, i) for i in label1]
fig_title1 = ["{}{}_tra_baselines.png".format(sample, i) for i in label1]

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.linewidth"] = 5
plt.rcParams["xtick.labelsize"] = 30
plt.rcParams["ytick.labelsize"] = 30
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["xtick.major.size"] = 12
plt.rcParams["ytick.major.size"] = 12
plt.rcParams["xtick.minor.size"] = 8
plt.rcParams["ytick.minor.size"] = 8
plt.rcParams["xtick.major.pad"] = 5
plt.rcParams["ytick.major.pad"] = 5
plt.rcParams["xtick.minor.pad"] = 3
plt.rcParams["ytick.minor.pad"] = 3
plt.rcParams["xtick.major.width"] = 5
plt.rcParams["ytick.major.width"] = 5
plt.rcParams["xtick.minor.width"] = 3
plt.rcParams["ytick.minor.width"] = 3

folder = "C:/python/calc_datasheets/{}/{}".format(date,sample)
pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
print(folder)

for run in range(0, len(label1)):

    df2 = pd.read_csv(datasheet2[run], low_memory=False)
    df_ini = pd.read_csv(datasheet2[run], low_memory=False)
    datalist2 = list(df2[:0])
    datalist_ext = datalist2[1:len(datalist2)-1]
    datalist_add = []
    for i in range(0,len(datalist_ext)):
        datalist_add_i = "{}_initial".format(datalist_ext[i])
        datalist_add.append(datalist_add_i)
    datalist_comb = datalist2 + datalist_add
    datalist_comb.insert(len(datalist2),"skip")
    datalist_comb.append("Time")
    # csv
    f1 = open("{}/{}".format(folder, csv_title1[run]), "w")
    writer1 = csv.writer(f1)
    csvlist = []
    for i in range(0, len(datalist2)-1):
        csvlist.append(datalist2[i])
    writer1.writerow(csvlist)
    f1.close()

    fig1 = plt.figure(figsize = (60, 10 * np.ceil(len(datalist2)/3)), linewidth=10, tight_layout=True)

    # search for outliers (NPs)
    DF1 = pd.read_csv("{}/{}".format(folder, csv_title1[run]))
    for i in range(0, len(datalist2)-1):
        DF1[datalist2[i]] = df2[datalist2[i]]

    # moving avg
    print("{}, bkg avg".format(sample))

    time = df2[datalist2[0]]
    iso_bkg_subt = []
    iso_bkg_subt_data = []
    total_NP_count = 0
    NP_count_hozon = pd.DataFrame()
    NP_count_hozon.insert(len(NP_count_hozon.columns),"",time)
    df_NP = pd.DataFrame()
    for i in range(1, len(datalist2)-1):
        isotope_i = DF1[datalist2[i]]
        isotope_i_data = df2[datalist2[i]]
        isotope_i_bkg = DF1[datalist2[i]]

        cycle = 0
        avg = []
        thres = []
        while True:
            m_avg = isotope_i_bkg.rolling(N, min_periods=1, center=True, win_type=None, on=None, axis=0, closed=None).mean()
            m_stdev = isotope_i_bkg.rolling(N, min_periods=1, center=True, win_type=None, on=None, axis=0, closed=None).std()
            thres_i = m_avg + n * m_stdev
            thres2_i = m_avg + n2 * m_stdev
            avg.append(m_avg)
            thres.append(thres_i)

            NP_discriminated = []
            bkg = []
            for k in range(0, len(isotope_i_bkg)):
                if isotope_i_bkg[k] <= thres_i[k]:
                    bkg.append(isotope_i_bkg[k])
                if isotope_i_bkg[k] > thres_i[k]:
                    NP_discriminated.append(isotope_i_bkg[k])
                    bkg.append(m_avg[k])
            isotope_i_bkg = pd.Series(bkg)
            cycle = cycle + 1

            if cycle == end:
                print("{}: {} counts".format(datalist2[i], np.round((sum(m_avg)/len(m_avg))*0.002, 2)))
                break

        # plot
        ax1 = fig1.add_subplot(np.ceil(len(datalist2)/3), 3, i)
        ax1.plot(time, isotope_i, marker="None", ls="-", lw=5, color="black", label="{}".format(datalist2[i].strip("'[]+")))

        cm = plt.cm.get_cmap("coolwarm")
        color_list = []
        A = len(avg)*2+1
        for idx in range(0,A):
            rgb = cm(idx/A)
            color_list.append(rgb)

        for X in range(0, len(avg)):
            ax1.plot(time, avg[X], color=color_list[len(avg)-X], ls="--", linewidth=10)
            ax1.plot(time, thres[X], color=color_list[len(avg)+1+X], ls="--", linewidth=10)

        ax1.set_xlabel("Time [msec]", fontsize=40, fontname="Arial", fontweight="bold")
        ax1.tick_params(axis="x", which="major", direction="out", length=15, pad=5, labelsize=40)
        ax1.set_ylabel("Signal Intensity [cps]", fontsize=40, fontname="Arial", fontweight="bold")
        ax1.tick_params(axis="y", which="major", direction="out", length=15, pad=5, labelsize=40)
        ax1.set_xlim(0, 30)
        ax1.set_ylim(0, thres[0].mean() * 3)
        ax1.grid()
        ax1.legend(fontsize=40, loc="upper right")
        ax1.set_axisbelow(True)

        NPs = 0
        for m in range(0, len(isotope_i)):
            if isotope_i[m] <= thres2_i[m] or (df2["TotalIonCurrent"][m] >= 3*10**7):
                isotope_i[m] = 0
                if isotope_i_data[m] <= thres2_i[m]:
                    isotope_i_data[m] = 0
                else:
                    isotope_i_data[m] = isotope_i_data[m] - m_avg[m]
            else:
                isotope_i[m] = isotope_i[m] - m_avg[m]
                isotope_i_data[m] = isotope_i_data[m] - m_avg[m]
                NPs = NPs + 1
        DF1[datalist2[i]] = isotope_i  # NP判別済 ＆ bkg減算済
        iso_bkg_subt.append(isotope_i)
        iso_bkg_subt_data.append(isotope_i_data)
        if datalist2[i] in isositei_list:
            total_NP_count = total_NP_count + isotope_i
            df_NP = pd.concat([df_NP,isotope_i],axis=1)
        NP_count_hozon.insert(len(NP_count_hozon.columns),i,isotope_i_data)

    DF1["TotalIonCurrent"] = df2["TotalIonCurrent"]

    NP_count_hozon.insert(len(NP_count_hozon.columns),"total",df2["TotalIonCurrent"])
    NP_count_hozon.columns = datalist2
    
    # fig output
    fig1.savefig("{}/{}".format(folder, fig_title1[run]))

    # csv output
    DF1.to_csv("{}/{}".format(folder, csv_title1[run]), index=None)

    NumofLines = mea_time*30040
    total_NP_count[NumofLines] = 0 

    skip = []
    skip_zengo = []
    skip_index = []
    NP_time_initial = []
    NP_time_interval = []
    for i in range(0, len(total_NP_count)):
        if i not in skip_index:
            add = 1
            if total_NP_count[i] > sum_thres:
                skip_i = [i]
                skip_zengo_i = [i]
                skip_index.append(i)
                NP_time_interval_i = [time[i]]
                NPiso_list = []
                for j in range(0,len(isositei_list)):
                    df_NP_j = df_NP[isositei_list[j]]
                    if not df_NP_j[i] == 0:
                        NPiso_list.append(isositei_list[j])
                while True:
                    if total_NP_count[i+add] > sum_thres:
                        NPiso_list_after = []
                        df_NP_j = df_NP[isositei_list[j]]
                        if not df_NP_j[i+add] == 0:
                            NPiso_list_after.append(isositei_list[j])
                        set_common = set(NPiso_list) & set(NPiso_list_after)
                        common_num = len(set_common)
                        if common_num == 0:
                            break
                        else:
                            skip_i.append(i+add)
                            skip_zengo_i.append(i+add)
                            skip_index.append(i+add)
                            NP_time_interval_i.append(time[i+add])
                            add = add + 1
                    if total_NP_count[i+add] <= sum_thres:
                        break
                skip.append(skip_i)
                skip_zengo.append(skip_zengo_i)
                NP_time_interval.append(NP_time_interval_i)
                NP_time_initial.append(NP_time_interval_i[0])

    print(len(df_ini)-5)

    for i in range(1,len(skip_zengo)):
        if skip_zengo[i][0] == 0:
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)
        elif skip_zengo[i][0] == 1:
            skip_first_1_i = int(skip[i][0])-1
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)
        elif skip_zengo[i][0] == 2:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)
        elif skip_zengo[i][0] == 3:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)
        elif skip_zengo[i][0] == 4:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)
        elif skip_zengo[i][len(skip[i])-1] == len(df_ini)-5:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
        elif skip_zengo[i][len(skip[i])-1] == len(df_ini)-4:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
        elif skip_zengo[i][len(skip[i])-1] == len(df_ini)-3:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
        elif skip_zengo[i][len(skip[i])-1] == len(df_ini)-2:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
            skip_zengo_i.append(skip_last_1_i)
        elif skip_zengo[i][len(skip[i])-1] == len(df_ini)-1:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
        else:
            skip_first_1_i = int(skip[i][0])-1
            skip_first_2_i = int(skip[i][0])-2
            skip_first_3_i = int(skip[i][0])-3
            skip_first_4_i = int(skip[i][0])-4
            skip_first_5_i = int(skip[i][0])-5
            skip_last_1_i = int(skip[i][len(skip[i])-1])+1
            skip_last_2_i = int(skip[i][len(skip[i])-1])+2
            skip_last_3_i = int(skip[i][len(skip[i])-1])+3
            skip_last_4_i = int(skip[i][len(skip[i])-1])+4
            skip_last_5_i = int(skip[i][len(skip[i])-1])+5
            skip_zengo_i = skip_zengo[i]
            skip_zengo_i.insert(0,skip_first_1_i)
            skip_zengo_i.insert(0,skip_first_2_i)
            skip_zengo_i.insert(0,skip_first_3_i)
            skip_zengo_i.insert(0,skip_first_4_i)
            skip_zengo_i.insert(0,skip_first_5_i)
            skip_zengo_i.append(skip_last_1_i)
            skip_zengo_i.append(skip_last_2_i)
            skip_zengo_i.append(skip_last_3_i)
            skip_zengo_i.append(skip_last_4_i)
            skip_zengo_i.append(skip_last_5_i)

    NP_iso = []
    total_iso = []
    NP_iso_data = []
    total_iso_data = []
    for i in range(0, len(iso_bkg_subt)):
        NP_iso_i = []
        total_iso_i = []
        NP_iso_data_i = []
        total_iso_data_i = []        
        for j in range(0, len(skip)):
            NP_iso_ij = []
            total_iso_ij = 0
            NP_iso_data_ij = []
            total_iso_data_ij = 0
            for k in skip[j]:
                NP_iso_ij.append(iso_bkg_subt[i][k])
                total_iso_ij = total_iso_ij + iso_bkg_subt[i][k]
                NP_iso_data_ij.append(iso_bkg_subt_data[i][k])
                total_iso_data_ij = total_iso_data_ij + iso_bkg_subt_data[i][k]
            NP_iso_i.append(NP_iso_ij)
            total_iso_i.append(total_iso_ij)
            NP_iso_data_i.append(NP_iso_data_ij)
            total_iso_data_i.append(total_iso_data_ij)
        NP_iso.append(NP_iso_i)
        total_iso.append(total_iso_i)
        NP_iso_data.append(NP_iso_data_i)
        total_iso_data.append(total_iso_data_i)

    csvlist_new = []
    for i in range(0, len(datalist2)-1):
        csvlist_new.append(datalist2[i])

    total_ion_curr = [df2["TotalIonCurrent"][i].values for i in skip]

    f5 = open("{}/{}".format(folder, csv_title5[run]), "w")
    writer5 = csv.writer(f5)
    csvlist_new_5 = []
    for i in range(0,len(datalist_comb)):
        csvlist_new_5.append(datalist_comb[i])
    writer5.writerow(csvlist_new_5)
    f5.close()
    
    DF5 = pd.read_csv("{}/{}".format(folder, csv_title5[run]))
    DF5[datalist_comb[0]] = pd.Series(NP_time_initial)

    for i in range(1, len(datalist2)-1):
        DF5[datalist2[i]] = pd.Series((np.array(total_iso_data[i-1])*int_time))

    total_ion_curr = [df2["TotalIonCurrent"][i].values for i in skip]
    DF5["TotalIonCurrent"] = total_ion_curr

    DF5["skip"] = skip_zengo

    for i in range(0,len(datalist_add)):
        data_initial_i = [df2[datalist2[i+1]][j].values for j in skip_zengo]
        DF5[datalist_add[i]] = data_initial_i

    time_all = [df2["t_elapsed_Buf"][i].values for i in skip_zengo]
    DF5["Time"] = time_all

    DF5.to_csv("{}/{}".format(folder, csv_title5[run]), index=None)

    print("done ({}/{})\n".format(run+1, len(label1)))
