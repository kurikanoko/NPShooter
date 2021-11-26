import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from numpy.core.arrayprint import dtype_is_implied
from numpy.lib.twodim_base import tri
import pandas as pd
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import csv
import pathlib
from tkinter import messagebox
import tkinter.simpledialog as simpledialog
import tkinter.filedialog
import pickle
import os
from tkinter import colorchooser
import sys
import glob
import re
from scipy.stats import f
from scipy import interpolate
import scipy.optimize
import math

#elementlist.binaryfileの場所を書き換える
f_binary=open("c:/python/elementlist.binaryfile","rb")
element_list = pickle.load(f_binary)

Atomic_Symbol_list = list(set(element_list["Atomic Symbol"])) 

#ソフトの色とフォント
color1="#ccccff" #背景
color2="#a3d1ff" #入力部分
color3="#b7ffb7" #三角グラフ背景
color4="#b7ffff" #三角柱グラフ背景
color5="#ffd6ad" #個数分布グラフ背景
color6="#ffbcff" #カウントグラフ背景
color7="#ffffc1" #新しいウィンドウ背景
color8="#ccffe5" #頂点編集するウィンドウ背景
label_font="Calibri"
label_font_jp="BIZ UDPゴシック"

#設定初期値
color_fig="white"
color_axes_pre="white"#グラフ表示していない状態での軸等の色
color_axes="black"#グラフ表示時の軸等の色
graph_font="Arial"#グラフのフォント
color_plot="blue"#プロットの色
alpha_plot=0.5 #プロットの透明度(0-1)
color_bar="green"#棒グラフの色
bar_style="枠あり"#棒グラフの見た目
weight_font="bold"
size_plot=35 #プロットのサイズ

plt.rcParams["font.weight"] = weight_font
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["axes.labelsize"] = 15
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5
plt.rcParams["xtick.major.pad"] = 2
plt.rcParams["ytick.major.pad"] = 2
plt.rcParams["xtick.major.width"] = 1.5
plt.rcParams["ytick.major.width"] = 1.5

iso_spe = ""

data_concat_syuturyoku_hanni = pd.DataFrame()
data_concat_syuturyoku_hanni_large = pd.DataFrame()

gosa_list_1 = []
gosa_list_2 = []
gosa_list_3 = []
gosa_list_4 = []
gosa_ox1 = ""
gosa_ox2 = ""
gosa_ox3 = ""
gosa_ox4 = ""

#参照値リスト
sansyou_1_a = ""
sansyou_1_b = ""
sansyou_1_c = ""
sansyou_1_sa = ""
sansyou_1_sb = ""
sansyou_1_sc = ""
sansyou_2_a = ""
sansyou_2_b = ""
sansyou_2_c = ""
sansyou_2_sa = ""
sansyou_2_sb = ""
sansyou_2_sc = ""
sansyou_3_a = ""
sansyou_3_b = ""
sansyou_3_c = ""
sansyou_3_sa = ""
sansyou_3_sb = ""
sansyou_3_sc = ""
sansyou_4_a = ""
sansyou_4_b = ""
sansyou_4_c = ""
sansyou_4_sa = ""
sansyou_4_sb = ""
sansyou_4_sc = ""

#範囲指定リスト
hanni_a_1 = ""
hanni_a_2 = ""
hanni_b_1 = ""
hanni_b_2 = ""
hanni_c_1 = ""
hanni_c_2 = ""
triangle_plot = ""

#頂点設定
#元素リスト
ele_p1 = "-"
ele_p2 = "-"
ele_p3 = "-"
ele_p4 = "-"
ele_p5 = "-"
ele_p6 = "-"
ele_p7 = "-"
ele_p8 = "-"
ele_p9 = "-"
ele_p10 = "-"
ele_p11 = "-"
ele_p12 = "-"
ele_p13 = "-"
ele_p14 = "-"
ele_p15 = "-"
ele_p_list = []
#質量数リスト
mass_p1 = "-"
mass_p2 = "-"
mass_p3 = "-"
mass_p4 = "-"
mass_p5 = "-"
mass_p6 = "-"
mass_p7 = "-"
mass_p8 = "-"
mass_p9 = "-"
mass_p10 = "-"
mass_p11 = "-"
mass_p12 = "-"
mass_p13 = "-"
mass_p14 = "-"
mass_p15 = "-"
mass_p_list = []
#補正項リスト
cor_p1 = ""
cor_p2 = ""
cor_p3 = ""
cor_p4 = ""
cor_p5 = ""
cor_p6 = ""
cor_p7 = ""
cor_p8 = ""
cor_p9 = ""
cor_p10 = ""
cor_p11 = ""
cor_p12 = ""
cor_p13 = ""
cor_p14 = ""
cor_p15 = ""
cor_p_list = []

#三次元組成分布色境
hi_bord_1 = 100
hi_bord_2 = 80
hi_bord_3 = 60
hi_bord_4 = 40
hi_bord_5 = 20

#誤差範囲パラメータ
i_number = 500

#三角柱グラフ縦粒径表示
size_hosei = ""
size_factor = ""
tate_iso = "all"
#三角柱グラフ縦粒径表示(頂点設定)
size_hosei_p = ""
size_factor_p = ""
tate_iso_p = "all"

#サイズ計算に使うかどうか(頂点設定)
size_calc_p1 = "true"
size_calc_p2 = "true"
size_calc_p3 = "true"
size_calc_p4 = "true"
size_calc_p5 = "true"
size_calc_p6 = "true"
size_calc_p7 = "true"
size_calc_p8 = "true"
size_calc_p9 = "true"
size_calc_p10 = "true"
size_calc_p11 = "true"
size_calc_p12 = "true"
size_calc_p13 = "true"
size_calc_p14 = "true"
size_calc_p15 = "true"

#三角グラフを着色するか
tri_plot_color = "OFF"

sample_label = ""

root=tk.Tk()
root.title(u"NP Shooter")
root.geometry("1850x980")
root.configure(bg=color1)

frame3=tk.Canvas(root,width=315,height=285,bg=color2)
frame3.place(x=8,y=28)

textframe=tk.Frame(root,width=150,height=35,bg=color2)
textframe.place(x=120,y=40)
textBox2=tk.Entry(textframe)
textBox2.place(width=145,height=20,x=0,y=0)
scrollbar=tk.Scrollbar(textframe,command=textBox2.xview,orient=tk.HORIZONTAL)
scrollbar.place(width=145,height=12,x=0,y=22)
textBox2["xscrollcommand"]=scrollbar.set

combobox2=ttk.Combobox(root,state="readonly",values=[])
combobox2.place(width=40,height=20,x=120,y=145)
combobox3=ttk.Combobox(root,state="readonly",values=[])
combobox3.place(width=40,height=20,x=170,y=145)
combobox4=ttk.Combobox(root,state="readonly",values=[])
combobox4.place(width=40,height=20,x=220,y=145)

combobox5=ttk.Combobox(root,state='readonly',values=[])
combobox5.place(width=40,height=20,x=120,y=100)
combobox6=ttk.Combobox(root,state='readonly',values=[])
combobox6.place(width=40,height=20,x=170,y=100)
combobox7=ttk.Combobox(root,state='readonly',values=[])
combobox7.place(width=40,height=20,x=220,y=100)

combobox8=ttk.Combobox(root,state='readonly',values=[])
combobox8.place(width=40,height=20,x=110,y=580)
combobox9=ttk.Combobox(root,state='readonly',values=[])
combobox9.place(width=40,height=20,x=110,y=610)


frame1=tk.Canvas(root,width=652,height=652,bg=color3)
frame1.place(x=334,y=22)
frame4=tk.Canvas(root,width=610,height=321,bg=color5)
frame4.place(x=996,y=22)
frame5=tk.Canvas(root,width=610,height=321,bg=color6)
frame5.place(x=996,y=353)
frame1_5=tk.Canvas(root,width=652,height=290,bg=color3)
frame1_5.place(x=334,y=684)

frame6=tk.Canvas(root,width=71,height=20,bg=color3)
frame6.place(x=16,y=326)
frame7=tk.Canvas(root,width=84,height=20,bg=color4)
frame7.place(x=16,y=436)
frame8=tk.Canvas(root,width=61,height=20,bg=color5)
frame8.place(x=16,y=546)
frame9=tk.Canvas(root,width=84,height=20,bg=color6)
frame9.place(x=16,y=686)

fig5=plt.figure(figsize=(6,3.11))
fig5.subplots_adjust(bottom=0.2,left=0.16)
ax5=fig5.add_subplot(111)
ax5.spines['top'].set_color("w")
ax5.spines['bottom'].set_color("w")
ax5.spines['left'].set_color("w")
ax5.spines['right'].set_color("w")
ax5.tick_params(colors="w")
canvas5=FigureCanvasTkAgg(fig5,master=root)
canvas5.get_tk_widget().place(x=1003,y=29)
canvas5._tkcanvas.place(x=1003,y=29)

fig6=plt.figure(figsize=(6,3.11))
fig6.subplots_adjust(bottom=0.2,left=0.16)
ax6=fig6.add_subplot(111)
ax6.spines['top'].set_color("w")
ax6.spines['bottom'].set_color("w")
ax6.spines['left'].set_color("w")
ax6.spines['right'].set_color("w")
ax6.tick_params(colors="w")
canvas6=FigureCanvasTkAgg(fig6,master=root)
canvas6.get_tk_widget().place(x=1003,y=360)
canvas6._tkcanvas.place(x=1003,y=360)

fig1_5=plt.figure(figsize=(6.42,2.8))
fig1_5.subplots_adjust(bottom=0.2,left=0.15,right=0.75,top=0.9)
ax1_5=fig1_5.add_subplot(111)
ax1_5.spines['top'].set_color("w")
ax1_5.spines['bottom'].set_color("w")
ax1_5.spines['left'].set_color("w")
ax1_5.spines['right'].set_color("w")
ax1_5.tick_params(colors="w")
canvas1_5=FigureCanvasTkAgg(fig1_5,master=root)
canvas1_5.get_tk_widget().place(x=341,y=691)
canvas1_5._tkcanvas.place(x=341,y=691)

def click_file():
    file_path = tkinter.filedialog.askopenfilename(filetypes = [("CSV(*.csv)","*.csv")])
    print(file_path)
    textBox2.delete(0, tkinter.END)
    textBox2.insert(tkinter.END,file_path)    


def click_element():
    global combobox5
    global combobox6   
    global combobox7
    global ele_list_min
    global ele_list
    global mass_list

    if not textBox2.get()=="":
        try:
            data_csv = pd.read_csv(textBox2.get())
            isotope_df_0=data_csv.columns
            data_csv_skip_column = int(data_csv.columns.get_loc("skip"))
            #isotope_df=isotope_df_0[1:len(isotope_df_0)-1]
            isotope_df=isotope_df_0[1:data_csv_skip_column-1]
            isotope_list=isotope_df.values.tolist()
            ele_list=[]
            for i in range(0,len(isotope_list)):
                alpha_i="".join([s for s in isotope_list[i] if s.isalpha()])
                ele_list.append(alpha_i)
            mass_list=[]
            for i in range(0,len(isotope_list)):
                digit_i="".join([s for s in isotope_list[i] if s.isdigit()])
                mass_list.append(digit_i)
            ele_list_min=list(set(ele_list))
            ele_list_min.sort()
            combobox5=ttk.Combobox(root,state='readonly',values=ele_list_min)
            combobox5.current(0)
            combobox5.place(width=40,height=20,x=120,y=100)
            combobox6=ttk.Combobox(root,state='readonly',values=ele_list_min)
            combobox6.current(0)
            combobox6.place(width=40,height=20,x=170,y=100)
            combobox7=ttk.Combobox(root,state='readonly',values=ele_list_min)
            combobox7.current(0)
            combobox7.place(width=40,height=20,x=220,y=100)
        except:
            messagebox.showerror('エラー','正しい形式のファイルを選択してください')    
    else:
        messagebox.showerror('エラー','ファイルを選択してください')

def click1():
   global combobox2
   global combobox3   
   global combobox4
   global mass_list_tB4
   global mass_list_tB5
   global mass_list_tB6

   if combobox5.get()==combobox6.get() or combobox6.get()==combobox7.get() or combobox5.get()==combobox7.get():
       messagebox.showerror('エラー','元素を重複なく選択してください')
   else:
       iso_list_tB4 = [k for k, x in enumerate(ele_list) if x == combobox5.get()]
       mass_list_tB4=[]
       for i in range(0,len(iso_list_tB4)):
           iso_index_i = iso_list_tB4[i]
           mass_i = mass_list[int(iso_index_i)]
           mass_list_tB4.append(mass_i)
       iso_list_tB5 = [k for k, x in enumerate(ele_list) if x == combobox6.get()]
       mass_list_tB5=[]
       for i in range(0,len(iso_list_tB5)):
           iso_index_i = iso_list_tB5[i]
           mass_i = mass_list[int(iso_index_i)]
           mass_list_tB5.append(mass_i)
       iso_list_tB6 = [k for k, x in enumerate(ele_list) if x == combobox7.get()]
       mass_list_tB6=[]
       for i in range(0,len(iso_list_tB6)):
           iso_index_i = iso_list_tB6[i]
           mass_i = mass_list[int(iso_index_i)]
           mass_list_tB6.append(mass_i)
       combobox2=ttk.Combobox(root,state='readonly',values=mass_list_tB4)
       combobox2.current(0)
       combobox2.place(width=40,height=20,x=120,y=145)
       combobox3=ttk.Combobox(root,state='readonly',values=mass_list_tB5)
       combobox3.current(0)
       combobox3.place(width=40,height=20,x=170,y=145)
       combobox4=ttk.Combobox(root,state='readonly',values=mass_list_tB6)
       combobox4.current(0)
       combobox4.place(width=40,height=20,x=220,y=145)

def click2():

    global combobox8
    global combobox9
    global ele_list_min

    if textBox2.get()=="":
        messagebox.showerror('エラー','ファイルを選択してください')

    if combobox2.get()=="":
        messagebox.showerror('エラー','質量数を選択してください')    
    else:
        label4_get=combobox2.get()
        print("label4:"+label4_get)

        label5_get=combobox3.get()
        print("label5:"+label5_get)

        label6_get=combobox4.get()
        print("label6:"+label6_get)

    sample_label_long = textBox2.get()
    global sample_label
    sample_label  = sample_label_long[:-22]

    all_file = glob.glob(sample_label+"_*_NP_events_large.csv")
    data_num = len(all_file)

    if textBox10.get() == "":
        count_min_X = 0
    else:
        count_min_X = textBox10.get()
    if textBox11.get() == "":
        count_min_Y = 0
    else:
        count_min_Y = textBox11.get()
    if textBox12.get() == "":
        count_min_Z = 0
    else:
        count_min_Z = textBox12.get()
    
    sample = ["{}_{}_NP_events_large".format(sample_label, str(i)) for i in range(1,data_num+1)]
    print("sample")
    print(sample)
    print(len(sample))

    global maximum_bar
    maximum_bar = len(sample)

    global data_large
    global data
    
    data_large = pd.DataFrame()
    for run in range(0, len(sample)):
        datasheet1 = "{}.csv".format(sample[run])
        print("\n{}".format(datasheet1))
        try:
            df1 = pd.read_csv(datasheet1, low_memory=False)
            data_large = data_large.append(df1, ignore_index=True)
        except FileNotFoundError:
            messagebox.showerror('エラー','選択されたファイルが正しくありません')

    data_skip_column = int(data_large.columns.get_loc("skip"))
    data = data_large.iloc[:,0:data_skip_column]
    global data_concat_syuturyoku
    data_concat_syuturyoku = data
    global data_concat_syuturyoku_large
    data_concat_syuturyoku_large = data_large

    tB4_data_1 = data["'[{}{}]+'".format(combobox2.get(),combobox5.get())]
    tB5_data_1 = data["'[{}{}]+'".format(combobox3.get(),combobox6.get())]
    tB6_data_1 = data["'[{}{}]+'".format(combobox4.get(),combobox7.get())]
    global data_gattai1
    data_gattai1=pd.concat([tB4_data_1,tB5_data_1,tB6_data_1],axis=1)
    data_gattai1.columns=["{}".format(combobox5.get()),"{}".format(combobox6.get()),"{}".format(combobox7.get())]
    gattai_drop_index = data_gattai1.index[(data_gattai1["{}".format(combobox5.get())] < float(count_min_X)) | (data_gattai1["{}".format(combobox6.get())] < float(count_min_Y)) | (data_gattai1["{}".format(combobox7.get())] < float(count_min_Z))]
    data_gattai1 = data_gattai1.drop(gattai_drop_index)

    data_concat_syuturyoku = data_concat_syuturyoku.drop(gattai_drop_index)
    data_concat_syuturyoku_large = data_concat_syuturyoku_large.drop(gattai_drop_index)

    global iso_compo_tB4
    global tB30
    global tB31
    global tB32
    iso_list_tB4 = [k for k, x in enumerate(element_list["Atomic Symbol"]) if x == combobox5.get()]
    mass_list_tB4=[]
    for i in range(0,len(iso_list_tB4)):
        iso_index_i = iso_list_tB4[i]
        massnumber = element_list["Mass Number"]
        mass_i = massnumber[int(iso_index_i)]
        mass_list_tB4.append(mass_i)
    iso_list_tB5 = [k for k, x in enumerate(element_list["Atomic Symbol"]) if x == combobox6.get()]
    mass_list_tB5=[]
    for i in range(0,len(iso_list_tB5)):
        iso_index_i = iso_list_tB5[i]
        massnumber = element_list["Mass Number"]
        mass_i = massnumber[int(iso_index_i)]
        mass_list_tB5.append(mass_i)
    iso_list_tB6 = [k for k, x in enumerate(element_list["Atomic Symbol"]) if x == combobox7.get()]
    mass_list_tB6=[]
    for i in range(0,len(iso_list_tB6)):
        iso_index_i = iso_list_tB6[i]
        massnumber = element_list["Mass Number"]
        mass_i = massnumber[int(iso_index_i)]
        mass_list_tB6.append(mass_i)

    Iso_Compo_list = element_list["Isotopic Composition"]
    list_number_tB4 = mass_list_tB4.index(int(combobox2.get()))
    index_tB4 = iso_list_tB4[int(list_number_tB4)]
    iso_compo_tB4 = Iso_Compo_list[int(index_tB4)]
    list_number_tB5 = mass_list_tB5.index(int(combobox3.get()))
    index_tB5 = iso_list_tB5[int(list_number_tB5)]
    iso_compo_tB5 = Iso_Compo_list[int(index_tB5)]
    list_number_tB6 = mass_list_tB6.index(int(combobox4.get()))
    index_tB6 = iso_list_tB6[int(list_number_tB6)]
    iso_compo_tB6 = Iso_Compo_list[int(index_tB6)]

    if textBox30.get()=="" and textBox31.get()=="" and textBox32.get()=="":
        textBox30.insert(tkinter.END,iso_compo_tB4)
        textBox31.insert(tkinter.END,iso_compo_tB5)
        textBox32.insert(tkinter.END,iso_compo_tB6)
        tB30=iso_compo_tB4
        tB31=iso_compo_tB5
        tB32=iso_compo_tB6
    else:
        tB30=textBox30.get()
        tB31=textBox31.get()
        tB32=textBox32.get()

    combobox8=ttk.Combobox(root,state='readonly',values=ele_list_min)
    combobox8.current(0)      
    combobox8.place(width=40,height=20,x=110,y=580)
    combobox9=ttk.Combobox(root,state='readonly',values=ele_list_min)
    combobox9.current(1)      
    combobox9.place(width=40,height=20,x=110,y=610)

    label11_2=tk.Label(root,text=combobox5.get()+":",bg=color1,font=(label_font,12))
    label11_2.place(width=20,height=20,x=80,y=720)
    label14_2=tk.Label(root,text=combobox6.get()+":",bg=color1,font=(label_font,12))
    label14_2.place(width=20,height=20,x=80,y=750)
    label17_2=tk.Label(root,text=combobox7.get()+":",bg=color1,font=(label_font,12))
    label17_2.place(width=20,height=20,x=80,y=780)

    ax3.cla()
    ax5.cla()
    ax6.cla()
    ax3.spines['top'].set_color(color_axes_pre)
    ax3.spines['bottom'].set_color(color_axes_pre)
    ax3.spines['left'].set_color(color_axes_pre)
    ax3.spines['right'].set_color(color_axes_pre)
    ax3.tick_params(colors=color_axes_pre)
    ax5.spines['top'].set_color(color_axes_pre)
    ax5.spines['bottom'].set_color(color_axes_pre)
    ax5.spines['left'].set_color(color_axes_pre)
    ax5.spines['right'].set_color(color_axes_pre)
    ax5.tick_params(colors=color_axes_pre)
    ax6.spines['top'].set_color(color_axes_pre)
    ax6.spines['bottom'].set_color(color_axes_pre)
    ax6.spines['left'].set_color(color_axes_pre)
    ax6.spines['right'].set_color(color_axes_pre)
    ax6.tick_params(colors=color_axes_pre)
    canvas3.draw()
    canvas5.draw()
    canvas6.draw()

    messagebox.showinfo('確認',"ファイル読み込みが完了しました")

#三角柱図
def clickA(): 

    plt.rcParams["font.family"] = graph_font
    
    #データ決め 
    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ,SUM],axis=1)
    data_concat.columns=['X','Y','Z','SUM']

    if size_hosei == "":
        data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
        data_z_x=data_gattai1[combobox5.get()]
    elif size_hosei == "sizefactor":
        if tate_iso == "all":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]*float(textBox30.get())/float(tB32))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        elif tate_iso == "A":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB30))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        else:
            messagebox.showerror("エラー","エラー")
    else:
        messagebox.showerror('エラー','エラー')

    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']

    h = np.sqrt(3.0)*0.5
   
    #座標・軸決め
    zmax=data_z['SUM'].max()

    amari=zmax%10

    if amari==0:
        baisuu=zmax
    else:
        baisuu=zmax-amari+10

    NP1=np.linspace(0,baisuu,6,dtype=int)

    plotx1=(100-data_concat['Y'])/100-data_concat['X']/200
    ploty1=h*data_concat['X']/100
    plotz1=data_z['SUM']/baisuu 
    data_plot1=pd.concat([plotx1,ploty1,plotz1],axis=1)
    data_plot1.columns=['x1','y1','z1']

    #プロット場所決め
    fig1=plt.figure(figsize=(8.5,7.5))
    ax1=fig1.add_subplot(111,projection='3d')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_zticks([])
    plt.axis('off')
    ax1.view_init(elev=12,azim=-69)

    #三角グラフ描き
    fig1.set_facecolor(color_fig)
    ax1.set_facecolor(color_fig)

    for i in range(1,5):
        ax1.plot([i*2/20.0, 1.0-i*2/20.0],[h*2*i/10.0, h*i*2/10.0],[0,0],color='gray', lw=0.5)
        ax1.plot([i*2/20.0, i*2/10.0],[h*i*2/10.0, 0.0],[0,0], color='gray', lw=0.5)
        ax1.plot([0.5+i*2/20.0, i*2/10.0],[h*(1.0-i*2/10.0), 0.0],[0,0],color='gray', lw=0.5)

    ax1.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
    ax1.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax1.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax1.plot([0.0, 0.0],[0.0, 0.0],[0,1], color=color_axes, lw=2)
    ax1.plot([1, 1],[0,0],[0,1], color=color_axes, lw=2)
    ax1.plot([0.5, 0.5],[h, h],[0,1], color=color_axes, lw=2)
    ax1.text(0.5, h+0.15,-0.1, combobox5.get(), fontsize=20,ha="center", color=color_axes)
    ax1.text(-0.15*h, -0.15/2,-0.1, combobox6.get(), fontsize=20,ha="center", color=color_axes)
    ax1.text(1+0.15*h, -0.15/2,-0.1, combobox7.get(), fontsize=20,ha="center", color=color_axes)

    for i in range(1,5):
        ax1.plot([2*i/20.0, 1.0-2*i/20.0],[2*h*i/10.0, 2*h*i/10.0],[1,1],color='gray', lw=0.5)
        ax1.plot([2*i/20.0, 2*i/10.0],[2*h*i/10.0, 0.0],[1,1], color='gray', lw=0.5)
        ax1.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[1,1], color='gray', lw=0.5)

    ax1.plot([0.0, 1.0],[0.0, 0.0],[1,1], color=color_axes, lw=2)
    ax1.plot([0.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)
    ax1.plot([1.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)

    ax1.plot([0,0],[0,0],[0,1], color=color_axes, lw=2)

    for i in range(0,6):
       ax1.plot([0,-0.02],[0,-0.02],[i/5,i/5], color=color_axes, lw=2)

    for i in range(0,6):
       ax1.text(-0.078,-0.078,i/5-0.01,NP1[i],fontsize=20,ha="right", color=color_axes)

    if size_hosei == "":
        ax1.text(-0.08,-0.08,1.2,"Counts",ha="right",fontsize=20,color=color_axes) 
    else:
        ax1.text(-0.08,-0.08,1.2,"Size (nm)",ha="right",fontsize=20,color=color_axes)
    
    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    ax1.set_zlim(0,1)

    #プロット
    ax1.scatter(data_plot1["x1"],data_plot1["y1"],data_plot1["z1"],c=color_plot,alpha=alpha_plot,depthshade=False,s=size_plot)

    #範囲指定のところに値が入っていた場合に全体表示のどの範囲にあたるかを示す
    if not textBox15.get() == "" or not textBox16.get() == "" :
        if textBox15.get()=="":
           tB15="0"
        else:
           tB15=textBox15.get()

        if textBox16.get()=="":
           tB16=baisuu
        else:
           tB16=textBox16.get()

        tB15_range = int(tB15)/baisuu
        tB16_range = int(tB16)/baisuu   

        ax1.plot([0.0, 1.0],[0.0, 0.0],[tB15_range,tB15_range], color="red", lw=1.5)
        ax1.plot([0.0, 0.5],[0.0, h],[tB15_range,tB15_range], color="red", lw=1.5)
        ax1.plot([1.0, 0.5],[0.0, h],[tB15_range,tB15_range], color="red", lw=1.5)
        ax1.plot([0.0, 0.0],[0.0, 0.0],[tB15_range,tB16_range], color="red", lw=1.5)
        ax1.plot([1, 1],[0,0],[tB15_range,tB16_range], color="red", lw=1.5)
        ax1.plot([0.5, 0.5],[h, h],[tB15_range,tB16_range], color="red", lw=1.5)
        ax1.plot([0.0, 1.0],[0.0, 0.0],[tB16_range,tB16_range], color="red", lw=1.5)
        ax1.plot([0.0, 0.5],[0.0, h],[tB16_range,tB16_range], color="red", lw=1.5)
        ax1.plot([1.0, 0.5],[0.0, h],[tB16_range,tB16_range], color="red", lw=1.5)        

    fig1.show()



#三角柱プロット範囲指定あり
def clickD(): 

    plt.rcParams["font.family"] = graph_font

    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ,SUM],axis=1)
    data_concat.columns=['X','Y','Z','SUM']

    if size_hosei == "":
        data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
        data_z_x=data_gattai1[combobox5.get()]
    elif size_hosei == "sizefactor":
        if tate_iso == "all":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]*float(textBox30.get())/float(tB31) +data_gattai1[combobox7.get()]*float(textBox30.get())/float(tB32))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        elif tate_iso == "A":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB30))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        else:
            messagebox.showerror("エラー","エラー")
    else:
        messagebox.showerror('エラー','エラー')

    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']

    h = np.sqrt(3.0)*0.5
     
    #座標決め
    zmax=data_z['SUM'].max()

    amari=zmax%10
    if amari==0:
       baisuu=zmax
    else:
       baisuu=zmax-amari+10

    plotx4=(100-data_concat['Y'])/100-data_concat['X']/200
    ploty4=h*data_concat['X']/100

    if textBox15.get()=="":
       tB15="0"
    else:
       tB15=textBox15.get()

    if textBox16.get()=="":
       tB16=baisuu
    else:
       tB16=textBox16.get()

    if tB15=="0" and tB16==baisuu:
       zvalue=baisuu
       plotz4=data_z['SUM']/zvalue
       NP4=np.linspace(0,baisuu,6,dtype=int)
    else:
       zvalue=int(tB16)-int(tB15)
       plotz4=(data_z['SUM']-int(tB15))/zvalue
       NP4=np.linspace(int(tB15),int(tB16),6,dtype=int)

    data_plot4=pd.concat([plotx4,ploty4,plotz4],axis=1)
    data_plot4.columns=['x4','y4','z4']
    data_plot40=data_plot4[(data_plot4["z4"]>=0)&(data_plot4["z4"]<=1)]

    #プロット場所決め
    fig4=plt.figure(figsize=(8.5,7.5))
    ax4=fig4.add_subplot(111,projection='3d')
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.set_zticks([])
    plt.axis('off')
    ax4.view_init(elev=12,azim=-69)

    #三角グラフ描き
    fig4.set_facecolor(color_fig)
    ax4.set_facecolor(color_fig)

    for i in range(1,5):
        ax4.plot([2*i/20.0, 1.0-2*i/20.0],[h*2*i/10.0, h*2*i/10.0],[0,0],color='gray', lw=0.5)
        ax4.plot([2*i/20.0, 2*i/10.0],[h*2*i/10.0, 0.0],[0,0], color='gray', lw=0.5)
        ax4.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[0,0], color='gray', lw=0.5)

    ax4.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
    ax4.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax4.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax4.plot([0.0, 0.0],[0.0, 0.0],[0,1], color=color_axes, lw=2)
    ax4.plot([1, 1],[0,0],[0,1], color=color_axes, lw=2)
    ax4.plot([0.5, 0.5],[h, h],[0,1], color=color_axes, lw=2)
    ax4.text(0.5, h+0.15,-0.1, combobox5.get(), fontsize=20,ha="center", color=color_axes)
    ax4.text(-0.15*h, -0.15/2,-0.1, combobox6.get(), fontsize=20,ha="center", color=color_axes)
    ax4.text(1+0.15*h, -0.15/2,-0.1, combobox7.get(), fontsize=20,ha="center", color=color_axes)

    for i in range(1,5):
        ax4.plot([2*i/20.0, 1.0-2*i/20.0],[2*h*i/10.0, 2*h*i/10.0],[1,1],color='gray', lw=0.5)
        ax4.plot([2*i/20.0, 2*i/10.0],[2*h*i/10.0, 0.0],[1,1], color='gray', lw=0.5)
        ax4.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[1,1], color='gray', lw=0.5)

    ax4.plot([0.0, 1.0],[0.0, 0.0],[1,1], color=color_axes, lw=2)
    ax4.plot([0.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)
    ax4.plot([1.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)

    ax4.plot([0,0],[0,0],[0,1], color=color_axes, lw=2)

    for i in range(0,6):
       ax4.plot([0,-0.02],[0,-0.02],[i/5,i/5], color=color_axes, lw=2)

    for i in range(0,6):
       ax4.text(-0.078,-0.078,i/5-0.01,NP4[i],fontsize=20,ha="right", color=color_axes)

    if size_hosei == "":
        ax4.text(-0.08,-0.08,1.2,"Counts",ha="right",fontsize=20,color=color_axes) 
    else:
        ax4.text(-0.08,-0.08,1.2,"Size (nm)",ha="right",fontsize=20,color=color_axes)
    ax4.set_xlim(0,1)
    ax4.set_ylim(0,1)
    ax4.set_zlim(0,1)

    #プロット
    ax4.scatter(data_plot40["x4"],data_plot40["y4"],data_plot40["z4"],c=color_plot,alpha=alpha_plot,depthshade=False,s=size_plot)

    fig4.show()


#三角プロット範囲指定あり
def clickC():

    plt.rcParams["font.family"] = graph_font

    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ,SUM],axis=1)
    data_concat.columns=['X','Y','Z','SUM']

    data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
    data_z_x=data_gattai1[combobox5.get()]
    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']           

    zmax=data_z['SUM'].max()

    print(data_z["SUM"].sort_values())

    amari=zmax%10
    if amari==0:
       baisuu=zmax
    else:
       baisuu=zmax-amari+10

    if textBox13.get()=="":
       tB13="0"
    else:
       tB13=textBox13.get()

    if textBox14.get()=="":
       tB14=baisuu
    else:
       tB14=textBox14.get()

    global data_size_click
    data_size_click = pd.concat([data_concat,data_z],axis=1)
    data_size_click.columns = ["a'","b'","c'","sum","SUM","X"]

    global triangle_plot
    triangle_plot = 1
    global data_concat0_tri
    if tB13=="0" and tB14==baisuu:
       data_concat0_tri=data_concat
       data_gattai0=data_gattai1
       data_size_click = data_size_click
       data_concat_syuturyoku_z = data_concat_syuturyoku
       data_concat_syuturyoku_z_large = data_concat_syuturyoku_large
    else:
       data_concat0_tri=data_concat[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]
       data_gattai0=data_gattai1[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]
       data_size_click = data_size_click[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]
       data_concat_syuturyoku_z = data_concat_syuturyoku[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]
       data_concat_syuturyoku_z_large = data_concat_syuturyoku_large[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]

    if textBox21.get()=="":
        tB21="0"
    else:
        tB21=textBox21.get()

    if textBox22.get()=="":
        tB22="100"
    else:
        tB22=textBox22.get()

    if textBox23.get()=="":
        tB23="0"
    else:
        tB23=textBox23.get()

    if textBox24.get()=="":
        tB24="100"
    else:
        tB24=textBox24.get()
        
    if textBox25.get()=="":
        tB25="0"
    else:
        tB25=textBox25.get()

    if textBox26.get()=="":
        tB26="100"
    else:
        tB26=textBox26.get()

    global data_concat_syuturyoku_hanni
    global data_concat_syuturyoku_hanni_large
    data_concat0_tri_hanni = data_concat0_tri[(data_concat0_tri["X"]>=int(tB21))&(data_concat0_tri["X"]<=int(tB22))&(data_concat0_tri["Y"]>=int(tB23))&(data_concat0_tri["Y"]<=int(tB24))&(data_concat0_tri["Z"]>=int(tB25))&(data_concat0_tri["Z"]<=int(tB26))]
    data_gattai_hanni = data_gattai0[(data_concat0_tri["X"]>=int(tB21))&(data_concat0_tri["X"]<=int(tB22))&(data_concat0_tri["Y"]>=int(tB23))&(data_concat0_tri["Y"]<=int(tB24))&(data_concat0_tri["Z"]>=int(tB25))&(data_concat0_tri["Z"]<=int(tB26))]
    data_size_click = data_size_click[(data_concat0_tri["X"]>=int(tB21))&(data_concat0_tri["X"]<=int(tB22))&(data_concat0_tri["Y"]>=int(tB23))&(data_concat0_tri["Y"]<=int(tB24))&(data_concat0_tri["Z"]>=int(tB25))&(data_concat0_tri["Z"]<=int(tB26))]
    data_concat_syuturyoku_hanni = data_concat_syuturyoku_z[(data_concat0_tri["X"]>=int(tB21))&(data_concat0_tri["X"]<=int(tB22))&(data_concat0_tri["Y"]>=int(tB23))&(data_concat0_tri["Y"]<=int(tB24))&(data_concat0_tri["Z"]>=int(tB25))&(data_concat0_tri["Z"]<=int(tB26))]
    data_concat_syuturyoku_hanni_large = data_concat_syuturyoku_z_large[(data_concat0_tri["X"]>=int(tB21))&(data_concat0_tri["X"]<=int(tB22))&(data_concat0_tri["Y"]>=int(tB23))&(data_concat0_tri["Y"]<=int(tB24))&(data_concat0_tri["Z"]>=int(tB25))&(data_concat0_tri["Z"]<=int(tB26))]

    #重み付けなし普通の平均とデータの標準誤差、標準誤差
    #平均値
    mean_data = data_concat0_tri_hanni.mean()
    print("平均")
    print(mean_data)
    #標準偏差
    std_data = np.std(data_concat0_tri_hanni, ddof=1) #標本標準偏差
    print("標準偏差")
    print(std_data)
    #標準誤差
    ste_data = std_data/(len(data_concat0_tri_hanni)**(1/2))
    print("標準誤差")
    print(ste_data)

    #重み付き平均値
    δ_data = data_gattai_hanni**(1/2)
    data_d_X = data_gattai_hanni[combobox5.get()]/float(textBox30.get())
    data_d_Y = data_gattai_hanni[combobox6.get()]/float(textBox31.get())
    data_d_Z = data_gattai_hanni[combobox7.get()]/float(textBox32.get())
    data_d_sum = data_d_X+data_d_Y+data_d_Z
    δ_data_d_X = δ_data[combobox5.get()]/float(textBox30.get())
    δ_data_d_Y = δ_data[combobox6.get()]/float(textBox31.get())
    δ_data_d_Z = δ_data[combobox7.get()]/float(textBox32.get())
    δ_s = (δ_data_d_X**2+δ_data_d_Y**2+δ_data_d_Z**2)**(1/2)
    δ_data_X = δ_data_d_X/data_d_X
    δ_data_X = δ_data_X.fillna(0)
    δ_data_Y = δ_data_d_Y/data_d_Y
    δ_data_Y = δ_data_Y.fillna(0)
    δ_data_Z = δ_data_d_Z/data_d_Z
    δ_data_Z = δ_data_Z.fillna(0)
    δ_X = 100*((δ_data_X)**2+(δ_s/data_d_sum)**2)**(1/2)
    δ_Y = 100*((δ_data_Y)**2+(δ_s/data_d_sum)**2)**(1/2)
    δ_Z = 100*((δ_data_Z)**2+(δ_s/data_d_sum)**2)**(1/2)
    σ_S = (δ_X**2+δ_Y**2+δ_Z**2)**(1/2)
    w_S = 1/((σ_S)**2)
    data_w_count_X = data_concat0_tri["X"]*w_S
    data_w_count_Y = data_concat0_tri["Y"]*w_S
    data_w_count_Z = data_concat0_tri["Z"]*w_S
    data_ave_w_X = np.sum(data_w_count_X)/np.sum(w_S)
    data_ave_w_Y = np.sum(data_w_count_Y)/np.sum(w_S)
    data_ave_w_Z = np.sum(data_w_count_Z)/np.sum(w_S)
    print("重みつき平均1")
    print("X")
    print(data_ave_w_X)
    print("Y")
    print(data_ave_w_Y)
    print("Z")
    print(data_ave_w_Z)
    #重み付き標準偏差
    data_zansa_X = data_concat0_tri["X"]-data_ave_w_X
    data_zansa_Y = data_concat0_tri["Y"]-data_ave_w_Y
    data_zansa_Z = data_concat0_tri["Z"]-data_ave_w_Z
    data_w_zansa2_X = w_S*(data_zansa_X**2)
    data_w_zansa2_Y = w_S*(data_zansa_Y**2)
    data_w_zansa2_Z = w_S*(data_zansa_Z**2)
    data_std_w_X = (np.sum(data_w_zansa2_X)/((len(data_concat)-1)*np.sum(w_S)))**(1/2)
    data_std_w_Y = (np.sum(data_w_zansa2_Y)/((len(data_concat)-1)*np.sum(w_S)))**(1/2)
    data_std_w_Z = (np.sum(data_w_zansa2_Z)/((len(data_concat)-1)*np.sum(w_S)))**(1/2)
    print("重みつき標準偏差1")
    print("X")
    print(data_std_w_X)
    print("Y")
    print(data_std_w_Y)
    print("Z")
    print(data_std_w_Z)

    h = np.sqrt(3.0)*0.5
    
    global data_plot3
    plotx3=(100-data_concat0_tri_hanni['Y'])/100-data_concat0_tri_hanni['X']/200
    ploty3=h*data_concat0_tri_hanni['X']/100
    data_plot3=pd.concat([plotx3,ploty3],axis=1)
    data_plot3.columns=['x3','y3']

    global data_click
    data_click=pd.concat([data_size_click,data_plot3],axis=1)
    data_click.columns=["a'","b'","c'","sum","SUM","X","x","y"]

    ax3.cla()
    fig3.set_facecolor(color_fig)
    ax3.set_facecolor(color_fig)

    for i in range(1,10):
        ax3.plot([i/20.0, 1.0-i/20.0],[h*i/10.0, h*i/10.0], linestyle='dashed',color='gray', lw=0.5)
        ax3.plot([i/20.0, i/10.0],[h*i/10.0, 0.0], linestyle='dashed',color='gray', lw=0.5)
        ax3.plot([0.5+i/20.0, i/10.0],[h*(1.0-i/10.0), 0.0], linestyle='dashed',color='gray', lw=0.5)

    ax3.plot([0.0, 1.0],[0.0, 0.0], color=color_axes, lw=2)
    ax3.plot([0.0, 0.5],[0.0, h], color=color_axes, lw=2)
    ax3.plot([1.0, 0.5],[0.0, h], color=color_axes, lw=2)
  
    ax3.text(0.455, h+0.0283, combobox5.get(), fontsize=22, color=color_axes)
    ax3.text(-0.1, -0.02, combobox6.get(), fontsize=22, color=color_axes)
    ax3.text(1.02, -0.02, combobox7.get(), fontsize=22, color=color_axes)

    for i in range(1,10):
        ax3.text(0.5+(10-i)/20.0+0.016, h*(1.0-(10-i)/10.0), '%d0' % i, fontsize=17, color=color_axes)
        ax3.text((10-i)/20.0-0.082, h*(10-i)/10.0, '%d0' % i, fontsize=17, color=color_axes)
        ax3.text(i/10.0-0.03, -0.06, '%d0' % i, fontsize=17, color=color_axes)

    ax3.text(-0.15,1,"Number of particles:"+str(len(data_plot3.dropna())),fontsize=14, color=color_axes)

    ax3.scatter(data_plot3["x3"],data_plot3["y3"],c=color_plot,alpha=alpha_plot,s=size_plot)

    if not sansyou_1_a =="":
        plot_sansyou_1_x = (100-float(sansyou_1_b))/100-float(sansyou_1_a)/200
        plot_sansyou_1_y = h*float(sansyou_1_a)/100
        global sansyou_1_sa
        global sansyou_1_sb
        global sansyou_1_sc
        if sansyou_1_sa == "":
            sansyou_1_sa = 0
        if sansyou_1_sb == "":
            sansyou_1_sb = 0
        if sansyou_1_sc == "":
            sansyou_1_sc = 0
        ax3.scatter(plot_sansyou_1_x,plot_sansyou_1_y,c="red",alpha=1,s=20)
        ax3.plot([plot_sansyou_1_x+float(sansyou_1_sa)*h/100/h/2,plot_sansyou_1_x-float(sansyou_1_sa)*h/100/h/2],[plot_sansyou_1_y-float(sansyou_1_sa)*h/100,plot_sansyou_1_y+float(sansyou_1_sa)*h/100],color="red",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_1_x-float(sansyou_1_sb)*h/100/h/2,plot_sansyou_1_x+float(sansyou_1_sb)*h/100/h/2],[plot_sansyou_1_y-float(sansyou_1_sb)*h/100,plot_sansyou_1_y+float(sansyou_1_sb)*h/100],color="red",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_1_x-float(sansyou_1_sc)*h/100/h,plot_sansyou_1_x+float(sansyou_1_sc)*h/100/h],[plot_sansyou_1_y,plot_sansyou_1_y],color="red",alpha=0.8,lw=1.2)

    if not sansyou_2_a =="":
        plot_sansyou_2_x = (100-float(sansyou_2_b))/100-float(sansyou_2_a)/200
        plot_sansyou_2_y = h*float(sansyou_2_a)/100
        global sansyou_2_sa
        global sansyou_2_sb
        global sansyou_2_sc
        if sansyou_2_sa == "":
            sansyou_2_sa = 0
        if sansyou_2_sb == "":
            sansyou_2_sb = 0
        if sansyou_2_sc == "":
            sansyou_2_sc = 0
        ax3.scatter(plot_sansyou_2_x,plot_sansyou_2_y,c="blue",alpha=1,s=20)
        ax3.plot([plot_sansyou_2_x+float(sansyou_2_sa)*h/100/h/2,plot_sansyou_2_x-float(sansyou_2_sa)*h/100/h/2],[plot_sansyou_2_y-float(sansyou_2_sa)*h/100,plot_sansyou_2_y+float(sansyou_2_sa)*h/100],color="blue",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_2_x-float(sansyou_2_sb)*h/100/h/2,plot_sansyou_2_x+float(sansyou_2_sb)*h/100/h/2],[plot_sansyou_2_y-float(sansyou_2_sb)*h/100,plot_sansyou_2_y+float(sansyou_2_sb)*h/100],color="blue",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_2_x-float(sansyou_2_sc)*h/100/h,plot_sansyou_2_x+float(sansyou_2_sc)*h/100/h],[plot_sansyou_2_y,plot_sansyou_2_y],color="blue",alpha=0.8,lw=1.2)

    if not sansyou_3_a =="":
        plot_sansyou_3_x = (100-float(sansyou_3_b))/100-float(sansyou_3_a)/200
        plot_sansyou_3_y = h*float(sansyou_3_a)/100
        global sansyou_3_sa
        global sansyou_3_sb
        global sansyou_3_sc
        if sansyou_3_sa == "":
            sansyou_3_sa = 0
        if sansyou_3_sb == "":
            sansyou_3_sb = 0
        if sansyou_3_sc == "":
            sansyou_3_sc = 0
        ax3.scatter(plot_sansyou_3_x,plot_sansyou_3_y,c="green",alpha=1,s=20)
        ax3.plot([plot_sansyou_3_x+float(sansyou_3_sa)*h/100/h/2,plot_sansyou_3_x-float(sansyou_3_sa)*h/100/h/2],[plot_sansyou_3_y-float(sansyou_3_sa)*h/100,plot_sansyou_3_y+float(sansyou_3_sa)*h/100],color="green",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_3_x-float(sansyou_3_sb)*h/100/h/2,plot_sansyou_3_x+float(sansyou_3_sb)*h/100/h/2],[plot_sansyou_3_y-float(sansyou_3_sb)*h/100,plot_sansyou_3_y+float(sansyou_3_sb)*h/100],color="green",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_3_x-float(sansyou_3_sc)*h/100/h,plot_sansyou_3_x+float(sansyou_3_sc)*h/100/h],[plot_sansyou_3_y,plot_sansyou_3_y],color="green",alpha=0.8,lw=1.2)
    
    if not sansyou_4_a =="":
        plot_sansyou_4_x = (100-float(sansyou_4_b))/100-float(sansyou_4_a)/200
        plot_sansyou_4_y = h*float(sansyou_4_a)/100
        global sansyou_4_sa
        global sansyou_4_sb
        global sansyou_4_sc
        if sansyou_4_sa == "":
            sansyou_4_sa = 0
        if sansyou_4_sb == "":
            sansyou_4_sb = 0
        if sansyou_4_sc == "":
            sansyou_4_sc = 0
        ax3.scatter(plot_sansyou_4_x,plot_sansyou_4_y,c="darkorange",alpha=1,s=20)
        ax3.plot([plot_sansyou_4_x+float(sansyou_4_sa)*h/100/h/2,plot_sansyou_4_x-float(sansyou_4_sa)*h/100/h/2],[plot_sansyou_4_y-float(sansyou_4_sa)*h/100,plot_sansyou_4_y+float(sansyou_4_sa)*h/100],color="darkorange",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_4_x-float(sansyou_4_sb)*h/100/h/2,plot_sansyou_4_x+float(sansyou_4_sb)*h/100/h/2],[plot_sansyou_4_y-float(sansyou_4_sb)*h/100,plot_sansyou_4_y+float(sansyou_4_sb)*h/100],color="darkorange",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_4_x-float(sansyou_4_sc)*h/100/h,plot_sansyou_4_x+float(sansyou_4_sc)*h/100/h],[plot_sansyou_4_y,plot_sansyou_4_y],color="darkorange",alpha=0.8,lw=1.2)

    canvas3.draw()

#組成分布
def clickE():
    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ,SUM],axis=1)
    data_concat.columns=['X','Y','Z','SUM']

    data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
    data_z_x=data_gattai1[combobox5.get()]
    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']           

    zmax=data_z['SUM'].max()

    amari=zmax%10
    if amari==0:
       baisuu=zmax
    else:
       baisuu=zmax-amari+10

    if textBox13.get()=="":
       tB13="0"
    else:
       tB13=textBox13.get()

    if textBox14.get()=="":
       tB14=baisuu
    else:
       tB14=textBox14.get()

    if tB13=="0" and tB14==baisuu:
       data_concat0=data_concat
    else:
       data_concat0=data_concat[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]


    if textBox18.get()=="":
        tB18="0"
    else:
        tB18=textBox18.get()

    if textBox19.get()=="":
        tB19="100"
    else:
        tB19=textBox19.get()

    if combobox8.get()==combobox5.get():
        if tB18=="0" and tB19=="100":
            data_concat00=data_concat0
        else:
            data_concat00=data_concat0[(data_concat0['X']>=int(tB18))&(data_concat0['X']<=int(tB19))]

        if combobox9.get()==combobox6.get():
            count=data_concat00['Y']/(data_concat00['Y']+data_concat00['Z'])*100
            perelement=combobox6.get()
        elif combobox9.get()==combobox7.get():
            count=data_concat00['Z']/(data_concat00['Y']+data_concat00['Z'])*100
            perelement=combobox7.get()   
        elif combobox9.get()=="":
            messagebox.showerror('エラー','割合を見る元素を指定してください') 
        else:
            messagebox.showerror('エラー','元素は'+combobox5.get()+"または"+combobox6.get()+"または"+combobox7.get()+'のいずれかを重複せずに入力してください') 
    elif combobox8.get()==combobox6.get():        
        if tB18=="0" and tB19=="100":
            data_concat00=data_concat0
        else:
            data_concat00=data_concat0[(data_concat0['Y']>=int(tB18))&(data_concat0['Y']<=int(tB19))]
      
        if combobox9.get()==combobox7.get():
            count=data_concat00['Z']/(data_concat00['Z']+data_concat00['X'])*100
            perelement=combobox7.get()
        elif combobox9.get()==combobox5.get():
            count=data_concat00['X']/(data_concat00['Z']+data_concat00['X'])*100
            perelement=combobox5.get() 
        elif combobox9.get()=="":
            messagebox.showerror('エラー','割合を見る元素を指定してください')       
        else:
            messagebox.showerror('エラー','元素は'+combobox5.get()+"または"+combobox6.get()+"または"+combobox7.get()+'のいずれかを入力してください') 
    elif combobox8.get()==combobox7.get():
        if tB18=="0" and tB19=="100":
            data_concat00=data_concat0
        else:
            data_concat00=data_concat0[(data_concat0['Z']>=int(tB18))&(data_concat0['Z']<=int(tB19))]
      
        if combobox9.get()==combobox5.get():
            count=data_concat00['X']/(data_concat00['X']+data_concat00['Y'])*100
            perelement=combobox5.get()
        elif combobox9.get()==combobox6.get():
            count=data_concat00['Y']/(data_concat00['X']+data_concat00['Y'])*100 
            perelement=combobox6.get()
        elif combobox9.get()=="":
            messagebox.showerror('error','割合を見る元素を指定してください') 
        else:
            messagebox.showerror('error','元素は'+combobox5.get()+"または"+combobox6.get()+"または"+combobox7.get()+'のいずれかを入力してください') 
    elif combobox8.get()=="":
        messagebox.showerror('error','範囲を固定する元素を指定してください')
    else:
        messagebox.showerror('error','元素は'+combobox5.get()+"または"+combobox6.get()+"または"+combobox7.get()+'のいずれかを入力してください')

    data_z_c=pd.concat([count,data_z_SUM],axis=1)
    data_z_c.columns=['X','SUM']     

    #重み付けなし普通の平均とデータの標準誤差、標準誤差
    #平均値
    mean_data = count.mean()
    print("平均")
    print(mean_data)
    #標準偏差
    std_data = np.std(count, ddof=1) #標本標準偏差
    print("標準偏差")
    print(std_data)
    #標準誤差
    ste_data = std_data/(len(count)**(1/2))
    print("標準誤差")
    print(ste_data)

    tB33=combobox1.get()

    times=100/int(tB33)

    NumX=[]
    for i in range(0,int(times)):
        NumX.append(np.count_nonzero((count>=i*int(tB33))&(count<=(i+1)*int(tB33))))

    ax5.cla()

    plt.rcParams["font.family"] = graph_font

    fig5.set_facecolor(color_fig)
    ax5.set_facecolor(color_fig)

    ax5.spines['top'].set_color(color_axes)
    ax5.spines['bottom'].set_color(color_axes)
    ax5.spines['left'].set_color(color_axes)
    ax5.spines['right'].set_color(color_axes)
    ax5.tick_params(colors=color_axes)

    left=np.linspace(0,100-int(tB33),int(times),dtype=int)
    height=np.array(NumX)
    if bar_style == "枠あり":   
        ax5.bar(left,height,width=int(tB33),color=color_bar,linewidth=1,edgecolor=color_axes,align="edge",zorder=2)
    elif bar_style == "枠なし":
        ax5.bar(left,height,width=int(tB33)* 0.8,color=color_bar,align="edge",zorder=2)
    ax5.grid(b=True,which='major',axis='y',color=color_axes,linewidth=0.5,zorder=1)
    ax5.set_xlabel("%"+perelement,fontname=graph_font,color=color_axes,weight=weight_font)
    ax5.set_ylabel("Number of particles",fontname=graph_font,color=color_axes,weight=weight_font)

    canvas5.draw()

#カウント分布
def clickF():

    plt.rcParams["font.family"] = graph_font

    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ],axis=1)
    data_concat.columns=['X','Y','Z']

    if size_hosei == "":
        data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
        data_z_x=data_gattai1[combobox5.get()]
    elif size_hosei == "sizefactor":
        if tate_iso == "all":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB31)+data_gattai1[combobox6.get()]*float(textBox30.get())/float(tB31) +data_gattai1[combobox7.get()]*float(textBox30.get())/float(tB32))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        elif tate_iso == "A":
            data_z_SUM=((data_gattai1[combobox5.get()]*float(textBox30.get())/float(tB30))*size_factor)**(1/3)
            data_z_x=data_gattai1[combobox5.get()]
        else:
            messagebox.showerror("エラー","エラー")
    else:
        messagebox.showerror('エラー','エラー')

    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']

    data_concat_z=pd.concat([data_concat,data_z],axis=1)

    if textBox21.get()=="":
        tB21="0"
    else:
        tB21=textBox21.get()

    if textBox22.get()=="":
        tB22="100"
    else:
        tB22=textBox22.get()

    if textBox23.get()=="":
        tB23="0"
    else:
        tB23=textBox23.get()

    if textBox24.get()=="":
        tB24="100"
    else:
        tB24=textBox24.get()
        
    if textBox25.get()=="":
        tB25="0"
    else:
        tB25=textBox25.get()

    if textBox26.get()=="":
        tB26="100"
    else:
        tB26=textBox26.get()
    
    data_concat0=data_concat_z[(data_concat_z["X"]>=int(tB21))&(data_concat_z["X"]<=int(tB22))&(data_concat_z["Y"]>=int(tB23))&(data_concat_z["Y"]<=int(tB24))&(data_concat_z["Z"]>=int(tB25))&(data_concat_z["Z"]<=int(tB26))]
    zmax0=data_concat0['SUM'].max()
    amari=zmax0%10
    if amari==0:
       baisuu=zmax0
    else:
       baisuu=zmax0-amari+10

    if textBox27.get()=="":
        tB27="0"
    else:
        tB27=textBox27.get()

    if textBox28.get()=="":
        tB28=baisuu
    else:
        tB28=textBox28.get()

    if textBox29.get()=="":
        tB29=100
        textBox29.insert(tkinter.END,100)
    else:
        tB29=textBox29.get()
    
    data_concat0_tri=data_concat0[(data_concat0["SUM"]>=int(tB27))&(data_concat0["SUM"]<=int(tB28))]

    if len(data_concat0_tri)==0:
        if int(tB25)<100-(int(tB22)+int(tB24)) or int(tB26)>100-(int(tB21)+int(tB23)):
            messagebox.showerror('エラー','範囲が間違っています')
        else:
            messagebox.showerror('エラー','範囲内にデータがありません')
    else:
        amari0_27=int(tB27)%int(tB29)
        if amari0_27==0:
            pltmin=int(tB27)
        else:
            pltmin=int(tB27)-amari0_27
        amari0_28=int(tB28)%int(tB29)
        if amari0_28==0:
            pltmax=int(tB28)
        else:
            pltmax=int(tB28)-amari0_28+int(tB29)

        ax6.cla()

        fig6.set_facecolor(color_fig)
        ax6.set_facecolor(color_fig)

        ax6.spines['top'].set_color(color_axes)
        ax6.spines['bottom'].set_color(color_axes)
        ax6.spines['left'].set_color(color_axes)
        ax6.spines['right'].set_color(color_axes)
        ax6.tick_params(colors=color_axes)
    
        if bar_style == "枠あり":
            ax6.hist(data_concat0_tri["SUM"],bins=np.arange(pltmin,pltmax+int(tB29),int(tB29)),color=color_bar,linewidth=1,edgecolor=color_axes,zorder=2)
        elif bar_style == "枠なし":
            ax6.hist(data_concat0_tri["SUM"],bins=np.arange(pltmin,pltmax+int(tB29),int(tB29)),width = int(tB29) * 0.8,color=color_bar,zorder=2)
        ax6.grid(b=True,which='major',axis='y',color=color_axes,linewidth=0.5,zorder=1)
        
        if size_hosei == "":
            ax6.set_xlabel("Counts",fontname=graph_font,color=color_axes,weight=weight_font)
        else:
            ax6.set_xlabel("Size (nm)",fontname=graph_font,color=color_axes,weight=weight_font)

        ax6.set_ylabel("Number of particles",fontname=graph_font,color=color_axes,weight=weight_font)

        data_concat0_small = data_concat0[(data_concat0["SUM"]>=int(tB27))&(data_concat0["SUM"]<=int(tB28))]

        ax6_pos = ax6.get_position()
        fig6.text(ax6_pos.x1-0.24,ax6_pos.y1+0.045,"{} particlespar".format(len(data_concat0)),fontsize=14, color=color_fig,backgroundcolor=color_fig)
        fig6.text(ax6_pos.x1-0.24,ax6_pos.y1+0.045,"{} particles".format(len(data_concat0_small)),fontsize=14, color=color_axes,backgroundcolor=color_fig)

        canvas6.draw()

#組成分布3D
def clickG():
    plt.rcParams["font.family"] = graph_font

    SUM=data_gattai1[combobox5.get()]/float(textBox30.get())+data_gattai1[combobox6.get()]/float(textBox31.get())+data_gattai1[combobox7.get()]/float(textBox32.get())
    pX=data_gattai1[combobox5.get()]/float(textBox30.get())/SUM*100 
    pY=data_gattai1[combobox6.get()]/float(textBox31.get())/SUM*100
    pZ=data_gattai1[combobox7.get()]/float(textBox32.get())/SUM*100
    data_concat=pd.concat([pX,pY,pZ,SUM],axis=1)
    data_concat.columns=['X','Y','Z','SUM']

    data_z_SUM=data_gattai1[combobox5.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB30)+data_gattai1[combobox6.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB31)+data_gattai1[combobox7.get()]/float(iso_compo_tB4)*float(textBox30.get())/float(tB32)
    data_z_x=data_gattai1[combobox5.get()] 
    data_z=pd.concat([data_z_SUM,data_z_x],axis=1)
    data_z.columns=['SUM','x']           

    zmax=data_z['SUM'].max()

    print(data_z["SUM"].sort_values())

    amari=zmax%10
    if amari==0:
       baisuu=zmax
    else:
       baisuu=zmax-amari+10

    if textBox13.get()=="":
       tB13="0"
    else:
       tB13=textBox13.get()

    if textBox14.get()=="":
       tB14=baisuu
    else:
       tB14=textBox14.get()

    global triangle_plot
    triangle_plot = 1
    global data_concat0_tri_g
    if tB13=="0" and tB14 == baisuu:
       data_concat0_tri_g = data_concat
    else:
       data_concat0_tri_g = data_concat[(data_z["SUM"]>=int(tB13))&(data_z["SUM"]<=int(tB14))]

    fig=plt.figure(figsize=(7.5,7.5))
    ax=fig.add_subplot(111,projection='3d')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.axis('off')
    ax.view_init(elev=12,azim=-69)

    h = np.sqrt(3.0)*0.5

    hi_list = []
    for i in range(0,10):
        for j in range(0,10-i):
          data_concat_hi_1 = data_concat0_tri_g[(data_concat0_tri_g["X"]>=i*10)&(data_concat0_tri_g["X"]<=(i+1)*10)&(data_concat0_tri_g["Y"]>=(9-i-j)*10)&(data_concat0_tri_g["Y"]<=(10-i-j)*10)&(data_concat0_tri_g["Z"]>=j*10)&(data_concat0_tri_g["Z"]<=(j+1)*10)]
          hi_1=len(data_concat_hi_1)
          hi_list.append(hi_1)

    for i in range(0,10):
        for j in range(0,10-i-1):
          data_concat_hi_2 = data_concat0_tri_g[(data_concat0_tri_g["X"]>=i*10)&(data_concat0_tri_g["X"]<=(i+1)*10)&(data_concat0_tri_g["Y"]>=(8-i-j)*10)&(data_concat0_tri_g["Y"]<=(9-i-j)*10)&(data_concat0_tri_g["Z"]>=j*10)&(data_concat0_tri_g["Z"]<=(j+1)*10)]
          hi_2=len(data_concat_hi_2)
          hi_list.append(hi_2)

    hi_max = max(hi_list)
    print("最大値")
    print(hi_max)

    for i in range(0,10):
        for j in range(0,10-i):
          data_concat_hi = data_concat0_tri_g[(data_concat0_tri_g["X"]>=i*10)&(data_concat0_tri_g["X"]<=(i+1)*10)&(data_concat0_tri_g["Y"]>=(9-i-j)*10)&(data_concat0_tri_g["Y"]<=(10-i-j)*10)&(data_concat0_tri_g["Z"]>=j*10)&(data_concat0_tri_g["Z"]<=(j+1)*10)]
          hi=len(data_concat_hi)
          if hi >= hi_bord_1:
              hi_color = "red"
          elif hi >= hi_bord_2:
              hi_color = "darkorange"
          elif hi >= hi_bord_3:
              hi_color = "lawngreen"
          elif hi >= hi_bord_4:
              hi_color = "green"
          elif hi >= hi_bord_5:
              hi_color = "blue"
          elif hi > 0:
              hi_color = "lightskyblue"
          else:
              hi_color = "#ffffff00"
          x1=[i/20+j/10,i/20+(j+1)/10,i/20+(j+1)/10,i/20+j/10]
          y1=[i*h/10,i*h/10,i*h/10,i*h/10]
          z1=[0,0,hi/hi_max,hi/hi_max]
          poly1=list(zip(x1,y1,z1))
          ax.add_collection3d(art3d.Poly3DCollection([poly1],facecolors=hi_color,linewidth=1,alpha=0.2))
          x2=[i/20+(j+1)/10,i/20+j/10+1/20,i/20+j/10+1/20,i/20+(j+1)/10]
          y2=[i*h/10,i*h/10+h/10,i*h/10+h/10,i*h/10]
          z2=[0,0,hi/hi_max,hi/hi_max]
          poly2=list(zip(x2,y2,z2))
          ax.add_collection3d(art3d.Poly3DCollection([poly2],color=hi_color,linewidth=1,alpha=0.2))
          x3=[i/20+j/10+1/20,i/20+j/10,i/20+j/10,i/20+j/10+1/20]
          y3=[i*h/10+h/10,i*h/10,i*h/10,i*h/10+h/10]
          z3=[0,0,hi/hi_max,hi/hi_max]
          poly3=list(zip(x3,y3,z3))
          ax.add_collection3d(art3d.Poly3DCollection([poly3],color=hi_color,linewidth=1,alpha=0.2))
          x4=[i/20+j/10,i/20+(j+1)/10,i/20+j/10+1/20]
          y4=[i*h/10,i*h/10,i*h/10+h/10]
          z4=[hi/hi_max,hi/hi_max,hi/hi_max]
          poly4=list(zip(x4,y4,z4))
          ax.add_collection3d(art3d.Poly3DCollection([poly4],color=hi_color,linewidth=1,alpha=0.2))
          x5=[i/20+j/10,i/20+(j+1)/10,i/20+j/10+1/20]
          y5=[i*h/10,i*h/10,i*h/10+h/10]
          z5=[0,0,0]
          poly5=list(zip(x5,y5,z5))
          ax.add_collection3d(art3d.Poly3DCollection([poly5],color=hi_color,linewidth=1,alpha=0.2))

    for i in range(0,10):
        for j in range(0,10-i-1):
          data_concat_hi = data_concat0_tri_g[(data_concat0_tri_g["X"]>=i*10)&(data_concat0_tri_g["X"]<=(i+1)*10)&(data_concat0_tri_g["Y"]>=(8-i-j)*10)&(data_concat0_tri_g["Y"]<=(9-i-j)*10)&(data_concat0_tri_g["Z"]>=j*10)&(data_concat0_tri_g["Z"]<=(j+1)*10)]
          hi=len(data_concat_hi)
          if hi >= hi_bord_1:
              hi_color = "red"
          elif hi >= hi_bord_2:
              hi_color = "darkorange"
          elif hi >= hi_bord_3:
              hi_color = "lawngreen"
          elif hi >= hi_bord_4:
              hi_color = "green"
          elif hi >= hi_bord_5:
              hi_color = "blue"
          elif hi > 0:
              hi_color = "lightskyblue"
          else:
              hi_color = "#ffffff00"
          x1=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+(j+1)/10,i/20+j/10+3/20]
          y1=[i*h/10+h/10,i*h/10,i*h/10,i*h/10+h/10]
          z1=[0,0,hi/hi_max,hi/hi_max]
          poly1=list(zip(x1,y1,z1))
          ax.add_collection3d(art3d.Poly3DCollection([poly1],facecolors=hi_color,linewidth=1,alpha=0.2))
          x2=[i/20+(j+1)/10,i/20+j/10+1/20,i/20+j/10+1/20,i/20+(j+1)/10]
          y2=[i*h/10,i*h/10+h/10,i*h/10+h/10,i*h/10]
          z2=[0,0,hi/hi_max,hi/hi_max]
          poly2=list(zip(x2,y2,z2))
          ax.add_collection3d(art3d.Poly3DCollection([poly2],color=hi_color,linewidth=1,alpha=0.2))
          x3=[i/20+j/10+1/20,i/20+j/10+3/20,i/20+j/10+3/20,i/20+j/10+1/20]
          y3=[i*h/10+h/10,i*h/10+h/10,i*h/10+h/10,i*h/10+h/10]
          z3=[0,0,hi/hi_max,hi/hi_max]
          poly3=list(zip(x3,y3,z3))
          ax.add_collection3d(art3d.Poly3DCollection([poly3],color=hi_color,linewidth=1,alpha=0.2))
          x4=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+j/10+1/20]
          y4=[i*h/10+h/10,i*h/10,i*h/10+h/10]
          z4=[hi/hi_max,hi/hi_max,hi/hi_max]
          poly4=list(zip(x4,y4,z4))
          ax.add_collection3d(art3d.Poly3DCollection([poly4],color=hi_color,linewidth=1,alpha=0.2))
          x5=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+j/10+1/20]
          y5=[i*h/10+h/10,i*h/10,i*h/10+h/10]
          z5=[0,0,0]
          poly5=list(zip(x5,y5,z5))
          ax.add_collection3d(art3d.Poly3DCollection([poly5],color=hi_color,linewidth=1,alpha=0.2))

    for i in range(1,10):
        ax.plot([2*i/40.0, 1.0-2*i/40.0],[h*2*i/20.0, h*2*i/20.0],[0,0],color='gray', lw=0.5)
        ax.plot([2*i/40.0, 2*i/20.0],[h*2*i/20.0, 0.0],[0,0], color='gray', lw=0.5)
        ax.plot([0.5+2*i/40.0, 2*i/20.0],[h*(1.0-2*i/20.0), 0.0],[0,0], color='gray', lw=0.5)

    ax.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
    ax.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
    ax.text(0.5, h+0.15,-0.1, combobox5.get(), fontsize=20,ha="center", color=color_axes)
    ax.text(-0.15*h, -0.15/2,-0.1, combobox6.get(), fontsize=20,ha="center", color=color_axes)
    ax.text(1+0.15*h, -0.15/2,-0.1, combobox7.get(), fontsize=20,ha="center", color=color_axes)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)

    fig.show()

def onclick(event):
    print("event.button=%d,event.xdata=%f,event.ydata=%f" %(event.button,event.xdata,event.ydata))
    eventx_min = event.xdata-0.02
    eventx_max = event.xdata+0.02
    eventy_min = event.ydata-0.02
    eventy_max = event.ydata+0.02
    data_click_kouho_1=data_click[(data_click["x"]>=eventx_min)&(data_click["x"]<=eventx_max)]
    data_click_kouho_2=data_click_kouho_1[(data_click_kouho_1["y"]>=eventy_min)&(data_click_kouho_1["y"]<=eventy_max)]
    data_y_sa = abs(data_click_kouho_2["y"]-event.ydata)
    data_click_satuki = pd.concat([data_click_kouho_2,data_y_sa],axis=1)
    data_click_satuki.columns = ["a'","b'","c'","sum","SUM","X","x","y","sa"]
    data_click_minjyun = data_click_satuki.sort_values("sa")
    data_click_minjyun_r = data_click_minjyun.reset_index()
    data_click_target_gyou = data_click_minjyun_r[0:1]
    print(data_click_target_gyou["index"])

    ax1_5.cla()

    plt.rcParams["font.family"] = graph_font

    fig1_5.set_facecolor(color_fig)
    ax1_5.set_facecolor(color_fig)

    ax1_5.spines['top'].set_color(color_axes)
    ax1_5.spines['bottom'].set_color(color_axes)
    ax1_5.spines['left'].set_color(color_axes)
    ax1_5.spines['right'].set_color(color_axes)
    ax1_5.tick_params(colors=color_axes)

    NP_number = int(data_click_target_gyou["index"])
    data_skip_column = int(data_large.columns.get_loc("skip"))
    element_number = int(data_skip_column)-2
    data_plot = data_large.iloc[:,int(data_skip_column):len(data_large)]
    data_time = data_large["Time"][NP_number]
    data_time_list = re.findall(r"\d+\.\d*",data_time)
    data_time_list_float = list(map(float,data_time_list))
    data_large_legend = data_large.columns[1:data_skip_column-1]

    def calc_002(n):
        return n*0.002

    if iso_spe == "":
        for i in range(0,element_number):
            data_element_i = data_large.iloc[NP_number,data_skip_column+i+1]
            data_element_list_i = re.findall(r"[+-]?[0-9]+\.[0-9]*[e]?[+-]?[0-9]*",data_element_i)
            data_element_list_float_ini_i = list(map(float,data_element_list_i))
            data_element_list_float_i = list(map(calc_002,data_element_list_float_ini_i))
            ax1_5.plot(data_time_list_float,data_element_list_float_i,lw=2)
    elif iso_spe == "specified":
        for i in range(0,element_number_spe):
            data_element_i = data_large.iloc[NP_number,data_skip_column+data_large_legend_spe_index[i]+1]
            data_element_list_i = re.findall(r"[+-]?[0-9]+\.[0-9]*[e]?[+-]?[0-9]*",data_element_i)
            data_element_list_float_ini_i = list(map(float,data_element_list_i))
            data_element_list_float_i = list(map(calc_002,data_element_list_float_ini_i))
            ax1_5.plot(data_time_list_float,data_element_list_float_i,lw=2)


    ax1_5.set_xlabel("Time(s)",fontname=graph_font,color=color_axes,weight=weight_font)
    ax1_5.set_ylabel("Intensity(Counts)",fontname=graph_font,color=color_axes,weight=weight_font)
    if iso_spe == "":
        ax1_5.legend(data_large_legend,bbox_to_anchor=(1.05,1),loc="upper left",borderaxespad=0)
    elif iso_spe == "specified":
        ax1_5.legend(data_large_legend_spe,bbox_to_anchor=(1.05,1),loc="upper left",borderaxespad=0)
    canvas1_5.draw()

    N = float(data_click_target_gyou["SUM"])

    p1_ = float(data_click_target_gyou["a'"])/100
    p2_ = float(data_click_target_gyou["b'"])/100
 
    X2 = 5.99

    def func_D(x):
        return (N*x*(x-2*p1_-2*p2_+2*p1_*p2_)+X2*x*(x-1)+N*p2_**2)**2-4*(-N*p2_**2+N*x+X2*x)*(N*p1_**2*x*(1-x))

    ti_x = (N*p2_**2)/(N+X2)

    y_p2_ = func_D(p2_)
    y_ti_x = func_D(ti_x)
    y_1 = func_D(1)

    if y_p2_*y_ti_x >= 0:
        left = scipy.optimize.fminbound(func_D,ti_x,p2_)
    else:
        left = ti_x

    if y_p2_*y_1 >= 0:
        right = scipy.optimize.fminbound(func_D,p2_,1)
    else:
        right = ti_x

    min_x = scipy.optimize.bisect(func_D,left,p2_)
    max_x = scipy.optimize.bisect(func_D,p2_,right)

    gosamaru_list = []
    for i in range(1,20):
        β_i = np.pi-(i*np.pi/20)
        p2_i = min_x+0.5*(np.cos(β_i)+1)*(max_x-min_x)
        a_i = -N*p2_**2+N*p2_i+X2*p2_i
        b_i = N*p2_i*(p2_i-2*p1_-2*p2_+2*p1_*p2_)+X2*p2_i*(p2_i-1)+N*p2_**2
        c_i = N*p1_**2*p2_i*(1-p2_i)
        D_i = b_i**2-4*a_i*c_i
        r_D_i = math.sqrt(D_i)
        j1_list_i = []
        j1_i = i
        j1_list_i.append(j1_i)
        min_p1_i = (-b_i-r_D_i)/(2*a_i)
        j1_list_i.append(min_p1_i)
        j1_list_i.append(p2_i)
        gosamaru_list.append(j1_list_i)
        j2_list_i = []
        j2_i = 40-i
        j2_list_i.append(j2_i)
        max_p1_i = (-b_i+math.sqrt(D_i))/(2*a_i)
        j2_list_i.append(max_p1_i)
        j2_list_i.append(p2_i)
        gosamaru_list.append(j2_list_i)


    a_20 = -N*p2_**2+N*min_x+X2*min_x
    b_20 = N*min_x*(min_x-2*p1_-2*p2_+2*p1_*p2_)+X2*min_x*(min_x-1)+N*p2_**2
    a_40 = -N*p2_**2+N*max_x+X2*max_x
    b_40 = N*max_x*(max_x-2*p1_-2*p2_+2*p1_*p2_)+X2*max_x*(max_x-1)+N*p2_**2
    j_list_20 = []
    j_20 = 40
    j_list_20.append(j_20)
    p1_20 = -b_20/(2*a_20)
    j_list_20.append(p1_20)
    j_list_20.append(min_x)
    gosamaru_list.append(j_list_20)
    j_list_40 = []
    j_40 = 20
    j_list_40.append(j_40)
    p1_40 = -b_40/(2*a_40)
    j_list_40.append(p1_40)
    j_list_40.append(max_x)
    gosamaru_list.append(j_list_40)

    β_1 = np.pi-(1*np.pi/20)
    p2_1 = min_x+0.5*(np.cos(β_1)+1)*(max_x-min_x)
    a_1 = -N*p2_**2+N*p2_1+X2*p2_1
    b_1 = N*p2_1*(p2_1-2*p1_-2*p2_+2*p1_*p2_)+X2*p2_1*(p2_1-1)+N*p2_**2
    c_1 = N*p1_**2*p2_1*(1-p2_1)
    D_1 = b_1**2-4*a_1*c_1
    r_D_1 = math.sqrt(D_1)
    j1_list_1 = []
    j1_1 = 41
    j1_list_1.append(j1_1)
    min_p1_1 = (-b_1-r_D_1)/(2*a_1)
    j1_list_1.append(min_p1_1)
    j1_list_1.append(p2_1)
    gosamaru_list.append(j1_list_1)
    gosamaru_df = pd.DataFrame(gosamaru_list,columns=["j","p1","p2"])
    gosamaru_df = gosamaru_df.sort_values("j")


    h = np.sqrt(3.0)*0.5

    gosamaru_plotx=(100-gosamaru_df["p2"]*100)/100-gosamaru_df["p1"]/2
    gosamaru_ploty=h*gosamaru_df["p1"]
    data_plot_gosa=pd.concat([gosamaru_plotx,gosamaru_ploty],axis=1)
    data_plot_gosa.columns=['p1','p2']

    ax3.cla()

    h=np.sqrt(3.0)*0.5

    for i in range(1,10):
        ax3.plot([i/20.0, 1.0-i/20.0],[h*i/10.0, h*i/10.0], linestyle='dashed',color='gray', lw=0.5)
        ax3.plot([i/20.0, i/10.0],[h*i/10.0, 0.0], linestyle='dashed',color='gray', lw=0.5)
        ax3.plot([0.5+i/20.0, i/10.0],[h*(1.0-i/10.0), 0.0], linestyle='dashed',color='gray', lw=0.5)

    ax3.plot([0.0, 1.0],[0.0, 0.0], color=color_axes, lw=2)
    ax3.plot([0.0, 0.5],[0.0, h], color=color_axes, lw=2)
    ax3.plot([1.0, 0.5],[0.0, h], color=color_axes, lw=2)
  
    ax3.text(0.455, h+0.0283, combobox5.get(), fontsize=22, color=color_axes)
    ax3.text(-0.1, -0.02, combobox6.get(), fontsize=22, color=color_axes)
    ax3.text(1.02, -0.02, combobox7.get(), fontsize=22, color=color_axes)

    for i in range(1,10):
        ax3.text(0.5+(10-i)/20.0+0.016, h*(1.0-(10-i)/10.0), '%d0' % i, fontsize=17, color=color_axes)
        ax3.text((10-i)/20.0-0.082, h*(10-i)/10.0, '%d0' % i, fontsize=17, color=color_axes)
        ax3.text(i/10.0-0.03, -0.06, '%d0' % i, fontsize=17, color=color_axes)

    ax3.text(-0.15,1,"Number of particles:"+str(len(data_plot3.dropna())),fontsize=14, color=color_axes)

    ax3.scatter(data_plot3["x3"],data_plot3["y3"],c=color_plot,alpha=alpha_plot,s=size_plot)

    if not sansyou_1_a =="":
        plot_sansyou_1_x = (100-float(sansyou_1_b))/100-float(sansyou_1_a)/200
        plot_sansyou_1_y = h*float(sansyou_1_a)/100
        global sansyou_1_sa
        global sansyou_1_sb
        global sansyou_1_sc
        if sansyou_1_sa == "":
            sansyou_1_sa = 0
        if sansyou_1_sb == "":
            sansyou_1_sb = 0
        if sansyou_1_sc == "":
            sansyou_1_sc = 0
        
        ax3.scatter(plot_sansyou_1_x,plot_sansyou_1_y,c="red",alpha=1,s=20)
        ax3.plot([plot_sansyou_1_x+float(sansyou_1_sa)*h/100/h/2,plot_sansyou_1_x-float(sansyou_1_sa)*h/100/h/2],[plot_sansyou_1_y-float(sansyou_1_sa)*h/100,plot_sansyou_1_y+float(sansyou_1_sa)*h/100],color="red",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_1_x-float(sansyou_1_sb)*h/100/h/2,plot_sansyou_1_x+float(sansyou_1_sb)*h/100/h/2],[plot_sansyou_1_y-float(sansyou_1_sb)*h/100,plot_sansyou_1_y+float(sansyou_1_sb)*h/100],color="red",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_1_x-float(sansyou_1_sc)*h/100/h,plot_sansyou_1_x+float(sansyou_1_sc)*h/100/h],[plot_sansyou_1_y,plot_sansyou_1_y],color="red",alpha=0.8,lw=1.2)

    if not sansyou_2_a =="":
        plot_sansyou_2_x = (100-float(sansyou_2_b))/100-float(sansyou_2_a)/200
        plot_sansyou_2_y = h*float(sansyou_2_a)/100
        global sansyou_2_sa
        global sansyou_2_sb
        global sansyou_2_sc
        if sansyou_2_sa == "":
            sansyou_2_sa = 0
        if sansyou_2_sb == "":
            sansyou_2_sb = 0
        if sansyou_2_sc == "":
            sansyou_2_sc = 0
        ax3.scatter(plot_sansyou_2_x,plot_sansyou_2_y,c="blue",alpha=1,s=20)
        ax3.plot([plot_sansyou_2_x+float(sansyou_2_sa)*h/100/h/2,plot_sansyou_2_x-float(sansyou_2_sa)*h/100/h/2],[plot_sansyou_2_y-float(sansyou_2_sa)*h/100,plot_sansyou_2_y+float(sansyou_2_sa)*h/100],color="blue",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_2_x-float(sansyou_2_sb)*h/100/h/2,plot_sansyou_2_x+float(sansyou_2_sb)*h/100/h/2],[plot_sansyou_2_y-float(sansyou_2_sb)*h/100,plot_sansyou_2_y+float(sansyou_2_sb)*h/100],color="blue",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_2_x-float(sansyou_2_sc)*h/100/h,plot_sansyou_2_x+float(sansyou_2_sc)*h/100/h],[plot_sansyou_2_y,plot_sansyou_2_y],color="blue",alpha=0.8,lw=1.2)

    if not sansyou_3_a =="":
        plot_sansyou_3_x = (100-float(sansyou_3_b))/100-float(sansyou_3_a)/200
        plot_sansyou_3_y = h*float(sansyou_3_a)/100
        global sansyou_3_sa
        global sansyou_3_sb
        global sansyou_3_sc
        if sansyou_3_sa == "":
            sansyou_3_sa = 0
        if sansyou_3_sb == "":
            sansyou_3_sb = 0
        if sansyou_3_sc == "":
            sansyou_3_sc = 0
        ax3.scatter(plot_sansyou_3_x,plot_sansyou_3_y,c="green",alpha=1,s=20)
        ax3.plot([plot_sansyou_3_x+float(sansyou_3_sa)*h/100/h/2,plot_sansyou_3_x-float(sansyou_3_sa)*h/100/h/2],[plot_sansyou_3_y-float(sansyou_3_sa)*h/100,plot_sansyou_3_y+float(sansyou_3_sa)*h/100],color="green",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_3_x-float(sansyou_3_sb)*h/100/h/2,plot_sansyou_3_x+float(sansyou_3_sb)*h/100/h/2],[plot_sansyou_3_y-float(sansyou_3_sb)*h/100,plot_sansyou_3_y+float(sansyou_3_sb)*h/100],color="green",alpha=0.8,lw=1.2)
        ax3.plot([plot_sansyou_3_x-float(sansyou_3_sc)*h/100/h,plot_sansyou_3_x+float(sansyou_3_sc)*h/100/h],[plot_sansyou_3_y,plot_sansyou_3_y],color="green",alpha=0.8,lw=1.2)

    ax3.scatter(data_click_target_gyou["x"],data_click_target_gyou["y"],c="orange",s=size_plot) 
    
    def spline(x,y,point,deg):
        tck,u = interpolate.splprep([x,y],k=deg,s=0) 
        u = np.linspace(0,1,num=point,endpoint=True) 
        spline = interpolate.splev(u,tck)
        return spline[0],spline[1]
    a,b = spline(data_plot_gosa["p1"],data_plot_gosa["p2"],100,3)
    ax3.plot(a,b,c="red")

    canvas3.draw()

fig3 = plt.figure(figsize=(6.42,6.42))
ax3 = fig3.add_subplot(111)
ax3.set_aspect('equal', 'datalim')
plt.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
plt.tick_params(bottom=False, left=False, right=False, top=False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
canvas3=FigureCanvasTkAgg(fig3,master=root)
canvas3.get_tk_widget().place(x=341,y=29)
canvas3._tkcanvas.place(x=341,y=29)
fig3.canvas.mpl_connect("button_press_event",onclick)

btn12=tk.Button(root,text="読み込み",command=click_element,font=(label_font_jp,8))
btn12.place(width=50,height=25,x=270,y=98)

btn4=tk.Button(root,text="表示",command=clickC,font=(label_font_jp,9))
btn4.place(width=50,height=25,x=215,y=359)

btn5=tk.Button(root,text="全体表示",command=clickA,font=(label_font_jp,9))
btn5.place(width=90,height=25,x=60,y=470)

btn6=tk.Button(root,text="表示",command=clickD,font=(label_font_jp,9))
btn6.place(width=50,height=25,x=215,y=499)

btn7=tk.Button(root,text="表示",command=clickE,font=(label_font_jp,9))
btn7.place(width=50,height=25,x=215,y=639)

btn8=tk.Button(root,text="表示",command=clickF,font=(label_font_jp,9))
btn8.place(width=50,height=25,x=215,y=839)

label24=tk.Label(root,text="分布",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label24.place(width=50,height=20,x=60,y=390)       

btn17=tk.Button(root,text="表示",command=clickG,font=(label_font_jp,9))
btn17.place(width=50,height=25,x=110,y=389)

btn2=tk.Button(root,text="csv読み込み",command=click2,font=(label_font_jp,9))
btn2.place(width=160,height=25,x=120,y=280)

label1=tk.Label(root,text="ファイル",anchor=tk.CENTER,bg=color2,font=(label_font_jp,9))
label1.place(width=100,height=20,x=12,y=45)
label3=tk.Label(root,text="元素",anchor=tk.CENTER,bg=color2,font=(label_font_jp,9))
label3.place(width=100,height=20,x=12,y=100)
label3_2=tk.Label(root,text="質量数",anchor=tk.CENTER,bg=color2,font=(label_font_jp,9))
label3_2.place(width=100,height=20,x=12,y=145)
label3_3=tk.Label(root,text="最低カウント数",anchor=tk.CENTER,bg=color2,font=(label_font_jp,9))
label3_3.place(width=100,height=20,x=12,y=190)
label3_4=tk.Label(root,text="補正",anchor=tk.CENTER,bg=color2,font=(label_font_jp,9))
label3_4.place(width=100,height=20,x=12,y=235)

btn11=tk.Button(root,text="参照",command=click_file,font=(label_font_jp,9))
btn11.place(width=40,height=25,x=275,y=43)
textBox10=tk.Entry(root)
textBox10.place(width=40,height=20,x=120,y=190)
textBox11=tk.Entry(root)
textBox11.place(width=40,height=20,x=170,y=190)
textBox12=tk.Entry(root)
textBox12.place(width=40,height=20,x=220,y=190)
textBox30=tk.Entry(root)
textBox30.place(width=40,height=20,x=120,y=235)
textBox31=tk.Entry(root)
textBox31.place(width=40,height=20,x=170,y=235)
textBox32=tk.Entry(root)
textBox32.place(width=40,height=20,x=220,y=235)

btn1=tk.Button(root,text="読み込み",command=click1,font=(label_font_jp,8))
btn1.place(width=50,height=25,x=270,y=143)

label20=tk.Label(root,text="三角グラフ",anchor=tk.W,bg=color3,font=(label_font_jp,10))
label20.place(width=67,height=16,x=20,y=330)

textBox13=tk.Entry(root)
textBox13.place(width=40,height=20,x=110,y=360)
textBox14=tk.Entry(root)
textBox14.place(width=40,height=20,x=160,y=360)
label7=tk.Label(root,text="-",bg=color1,font=(label_font,10))
label7.place(width=10,height=20,x=150,y=360)
label24=tk.Label(root,text="カウント",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label24.place(width=50,height=20,x=60,y=360)       

label21=tk.Label(root,text="三角柱グラフ",anchor=tk.W,bg=color4,font=(label_font_jp,10))
label21.place(width=80,height=16,x=20,y=440)

textBox15=tk.Entry(root)
textBox15.place(width=40,height=20,x=110,y=500)
textBox16=tk.Entry(root)
textBox16.place(width=40,height=20,x=160,y=500)
label8=tk.Label(root,text="-",bg=color1,font=(label_font,10))
label8.place(width=10,height=20,x=150,y=500)
label25=tk.Label(root,text="カウント",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label25.place(width=50,height=20,x=60,y=500)       

label22=tk.Label(root,text="組成分布",anchor=tk.W,bg=color5,font=(label_font_jp,10))
label22.place(width=57,height=16,x=20,y=550)

textBox18=tk.Entry(root)
textBox18.place(width=40,height=20,x=156,y=580)
textBox19=tk.Entry(root)
textBox19.place(width=40,height=20,x=205,y=580)
combobox1=ttk.Combobox(root,state='readonly',values=[1,2,4,5,10,20,25,50,100])
combobox1.current(4)      
combobox1.place(width=40,height=20,x=110,y=640)
label9=tk.Label(root,text=":",bg=color1,font=(label_font_jp,10))
label9.place(width=3,height=20,x=151,y=580)
label10=tk.Label(root,text="-",bg=color1,font=(label_font_jp,10))
label10.place(width=10,height=20,x=195,y=580)
label29=tk.Label(root,text="%",bg=color1,font=(label_font_jp,9))
label29.place(width=20,height=20,x=248,y=580) 
label26=tk.Label(root,text="固定",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label26.place(width=50,height=20,x=60,y=580) 
label27=tk.Label(root,text="表示",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label27.place(width=50,height=20,x=60,y=610) 
label30=tk.Label(root,text="幅",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label30.place(width=50,height=20,x=60,y=640) 

label22=tk.Label(root,text="カウント分布",anchor=tk.W,bg=color6,font=(label_font_jp,10))
label22.place(width=80,height=16,x=20,y=690)

textBox21=tk.Entry(root)
textBox21.place(width=40,height=20,x=110,y=720)
textBox22=tk.Entry(root)
textBox22.place(width=40,height=20,x=160,y=720)
textBox23=tk.Entry(root)
textBox23.place(width=40,height=20,x=110,y=750)
textBox24=tk.Entry(root)
textBox24.place(width=40,height=20,x=160,y=750)
textBox25=tk.Entry(root)
textBox25.place(width=40,height=20,x=110,y=780)
textBox26=tk.Entry(root)
textBox26.place(width=40,height=20,x=160,y=780)
label11=tk.Label(root,text=combobox5.get()+":",bg=color1,font=(label_font_jp,10))
label11.place(width=20,height=20,x=80,y=720)
label12=tk.Label(root,text="-",bg=color1,font=(label_font_jp,10))
label12.place(width=10,height=20,x=150,y=720)
label13=tk.Label(root,text="%",bg=color1,font=(label_font_jp,9))
label13.place(width=20,height=20,x=203,y=720)
label14=tk.Label(root,text=combobox6.get()+":",bg=color1,font=(label_font_jp,10))
label14.place(width=20,height=20,x=80,y=750)
label15=tk.Label(root,text="-",bg=color1,font=(label_font_jp,10))
label15.place(width=10,height=20,x=150,y=750)
label16=tk.Label(root,text="%",bg=color1,font=(label_font_jp,9))
label16.place(width=20,height=20,x=203,y=750)
label17=tk.Label(root,text=combobox7.get()+":",bg=color1,font=(label_font_jp,10))
label17.place(width=20,height=20,x=80,y=780)
label18=tk.Label(root,text="-",bg=color1,font=(label_font_jp,10))
label18.place(width=10,height=20,x=150,y=780)
label19=tk.Label(root,text="%",bg=color1,font=(label_font_jp,9))
label19.place(width=20,height=20,x=203,y=780)

textBox27=tk.Entry(root)
textBox27.place(width=40,height=20,x=110,y=810)
textBox28=tk.Entry(root)
textBox28.place(width=40,height=20,x=160,y=810)      
label28=tk.Label(root,text="-",bg=color1,font=(label_font,10))
label28.place(width=10,height=20,x=150,y=810)
label31=tk.Label(root,text="カウント",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label31.place(width=50,height=20,x=60,y=810) 
textBox29=tk.Entry(root)
textBox29.place(width=40,height=20,x=110,y=840) 
label32=tk.Label(root,text="幅",anchor=tk.CENTER,bg=color1,font=(label_font_jp,9))
label32.place(width=50,height=20,x=60,y=840) 

def savelist():
    tBlist = ["",textBox2.get(),"",combobox5.get(),combobox6.get(),combobox7.get(),combobox2.get(),combobox3.get(),combobox4.get(),textBox10.get(),textBox11.get(),textBox12.get(),textBox13.get(),textBox14.get(),textBox15.get(),textBox16.get(),combobox8.get(),textBox18.get(),textBox19.get(),combobox9.get(),textBox21.get(),textBox22.get(),textBox23.get(),textBox24.get(),textBox25.get(),textBox26.get(),textBox28.get(),textBox29.get(),textBox30.get(),textBox31.get(),textBox32.get(),combobox1.get()]
    tBlist_str = "?".join(tBlist)

    file_path = tkinter.filedialog.asksaveasfilename(defaultextension="txt",filetypes = [("TRIファイル(*.tri)","*.tri")])

    if len(file_path) != 0:
        f = open(file_path,"w")
        f.write(tBlist_str)
        f.close
    
def openlist():
    file_path = tkinter.filedialog.askopenfilename(filetypes = [("TRIファイル(*.tri)","*.tri")])

    if len(file_path) != 0:
        f = open(file_path)
        tBlist_txt = f.read()
        f.close()

        str = tBlist_txt
        tBlist = str.split("?")

        global textBox2
        textBox2.delete(0, tkinter.END)
        textBox2.insert(tkinter.END,tBlist[1])
        click_element()

        global combobox5
        global combobox6   
        global combobox7

        combo5_list_number=ele_list_min.index(tBlist[3])
        combobox5.current(int(combo5_list_number))

        combo6_list_number=ele_list_min.index(tBlist[4])
        combobox6.current(int(combo6_list_number))

        combo7_list_number=ele_list_min.index(tBlist[5])
        combobox7.current(int(combo7_list_number))

        click1()

        global combobox2
        global combobox3   
        global combobox4

        global combobox8
        global combobox9

        combo2_list_number=mass_list_tB4.index(tBlist[6])
        combobox2.current(int(combo2_list_number))

        combo3_list_number=mass_list_tB5.index(tBlist[7])
        combobox3.current(int(combo3_list_number))

        combo4_list_number=mass_list_tB6.index(tBlist[8])
        combobox4.current(int(combo4_list_number))

        textBox10.delete(0, tkinter.END)
        textBox10.insert(tkinter.END,tBlist[9])
        textBox11.delete(0, tkinter.END)
        textBox11.insert(tkinter.END,tBlist[10])
        textBox12.delete(0, tkinter.END)
        textBox12.insert(tkinter.END,tBlist[11])
        textBox13.delete(0, tkinter.END)
        textBox13.insert(tkinter.END,tBlist[12])
        textBox14.delete(0, tkinter.END)
        textBox14.insert(tkinter.END,tBlist[13])
        textBox15.delete(0, tkinter.END)
        textBox15.insert(tkinter.END,tBlist[14])
        textBox16.delete(0, tkinter.END)
        textBox16.insert(tkinter.END,tBlist[15])
        textBox18.delete(0, tkinter.END)
        textBox18.insert(tkinter.END,tBlist[17])
        textBox19.delete(0, tkinter.END)
        textBox19.insert(tkinter.END,tBlist[18])
        textBox21.delete(0, tkinter.END)
        textBox21.insert(tkinter.END,tBlist[20])
        textBox22.delete(0, tkinter.END)
        textBox22.insert(tkinter.END,tBlist[21])
        textBox23.delete(0, tkinter.END)
        textBox23.insert(tkinter.END,tBlist[22])
        textBox24.delete(0, tkinter.END)
        textBox24.insert(tkinter.END,tBlist[23])
        textBox25.delete(0, tkinter.END)
        textBox25.insert(tkinter.END,tBlist[24])
        textBox26.delete(0, tkinter.END)
        textBox26.insert(tkinter.END,tBlist[25])
        textBox28.delete(0, tkinter.END)
        textBox28.insert(tkinter.END,tBlist[26])
        textBox29.delete(0, tkinter.END)
        textBox29.insert(tkinter.END,tBlist[27])
        textBox30.delete(0, tkinter.END)
        textBox30.insert(tkinter.END,tBlist[28])
        textBox31.delete(0, tkinter.END)
        textBox31.insert(tkinter.END,tBlist[29])
        textBox32.delete(0, tkinter.END)
        textBox32.insert(tkinter.END,tBlist[30])
        combo1_list=[1,2,4,5,10,20,25,50,100]
        combo1_list_number=combo1_list.index(int(tBlist[31]))
        combobox1.current(int(combo1_list_number))

        click2()
        combo8_list_number=ele_list_min.index(tBlist[16])
        combobox8.current(int(combo8_list_number))
        combo9_list_number=ele_list_min.index(tBlist[19])
        combobox9.current(int(combo9_list_number))

        clickC()
        clickE()
        clickF()

def savefig3():
    file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
    
    print(file_path)

    if len(file_path) != 0:
        fig3.savefig(file_path)


def savefig5():
    file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
    
    print(file_path)

    if len(file_path) != 0:
        fig5.savefig(file_path)

def savefig6():
    file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
    
    print(file_path)

    if len(file_path) != 0:
        fig6.savefig(file_path)

def save_csv_large():
    if sample_label == "":
        messagebox.showerror("エラー","三角グラフを表示してください")
    else:
        print(sample_label)
        newwindow_large = tk.Toplevel(root)
        newwindow_large.attributes("-topmost",True)
        newwindow_large.geometry("500x120")
        newwindow_large.title(u"範囲指定csv出力(large)")
        label_large_1 = tk.Label(newwindow_large,text="ファイル保存",font=(label_font_jp,10),anchor=tk.W)
        label_large_1.place(width=80,height=20,x=10,y=10)
        label_large_2 = tk.Label(newwindow_large,text=sample_label,font=(label_font_jp,10),anchor=tk.W)
        label_large_2.place(width=480,height=20,x=20,y=30)
        label_large_3 = tk.Label(newwindow_large,text="_",font=(label_font_jp,10),anchor=tk.W)
        label_large_3.place(width=10,height=20,x=20,y=50)
        textBox_large_1 = ttk.Entry(newwindow_large)
        textBox_large_1.place(width=100,height=20,x=28,y=50)
        label_large_4 = tk.Label(newwindow_large,text="_1_NP_events_large.csv",font=(label_font_jp,10),anchor=tk.W)
        label_large_4.place(width=200,height=20,x=130,y=50)

        def OK_large():
            filename_large = textBox_large_1.get()
            data_concat_syuturyoku_hanni_large.to_csv("{}_{}_1_NP_events_large.csv".format(sample_label,filename_large))
            newwindow_large.destroy()

        def Cancel_large():
            newwindow_large.destroy()

        btn_large_1=tk.Button(newwindow_large,text="決定",command=OK_large)
        btn_large_1.place(width=80,height=25,x=310,y=80)
        btn_large_2=tk.Button(newwindow_large,text="キャンセル",command=Cancel_large)
        btn_large_2.place(width=80,height=25,x=400,y=80)

def layout():
    newwindow = tk.Toplevel(root)
    newwindow.geometry("500x500")
    newwindow.title(u"レイアウト設定")

    label33=tk.Label(newwindow,text="グラフのフォント",font=(label_font_jp,9),anchor=tk.W)
    label33.place(width=100,height=20,x=20,y=80)
    label34=tk.Label(newwindow,text="プロットの大きさ",font=(label_font_jp,9),anchor=tk.W)
    label34.place(width=100,height=20,x=20,y=110)
    label35=tk.Label(newwindow,text="プロットの透明度",font=(label_font_jp,9),anchor=tk.W)
    label35.place(width=100,height=20,x=20,y=140)
    label36=tk.Label(newwindow,text="プロットの色",font=(label_font_jp,9),anchor=tk.W)
    label36.place(width=100,height=20,x=20,y=170)
    label37=tk.Label(newwindow,text="棒グラフの色",font=(label_font_jp,9),anchor=tk.W)
    label37.place(width=100,height=20,x=20,y=200)

    bln1 = tkinter.BooleanVar()
    if color_fig == "black":
        bln1.set(True)
    elif color_fig == "white":
        bln1.set(False)
    else:
        print("checkbox_error")   
    check1 = tk.Checkbutton(newwindow,variable=bln1,text="背景を黒にする",font=(label_font_jp,9))
    check1.place(x=20,y=20)

    graph_font_list=["Arial","Calibri","Cambria","Meiryo"]
    graph_font_list_number=graph_font_list.index(graph_font)
    combobox12=ttk.Combobox(newwindow,state="readonly",values=["Arial","Calibri","Cambria","Meiryo"])
    combobox12.place(width=80,height=20,x=130,y=80)
    combobox12.current(int(graph_font_list_number))    

    color_plot_list_0=["red","blue","green","yellow","gray"]
    if color_plot in color_plot_list_0:
        color_plot_list = ["red","blue","green","yellow","gray"]
    else:
        color_plot_list = ["red","blue","green","yellow","gray",color_plot]
    color_plot_list_number=color_plot_list.index(color_plot)
    global combobox13
    combobox13=ttk.Combobox(newwindow,state="readonly",values=color_plot_list)
    combobox13.place(width=80,height=20,x=130,y=170)
    combobox13.current(int(color_plot_list_number))  

    def color_choose_plot():
        color_now = color_plot
        color = colorchooser.askcolor(color_now,title="その他の色",master=newwindow)
        color_list=list(color)
        colorcode=color_list[1]
        if not color == (None, None):
            global combobox13
            combobox13=ttk.Combobox(newwindow,state="readonly",values=["red","blue","green","yellow","gray",colorcode])
            combobox13.place(width=80,height=20,x=130,y=170)
            combobox13.current(5)  
        newwindow.attributes("-topmost",True)            

    btn11=tk.Button(newwindow,text="その他の色",command=color_choose_plot)
    btn11.place(width=80,height=25,x=225,y=170)    

    alpha_plot_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    alpha_plot_list_number=alpha_plot_list.index(alpha_plot)
    combobox14=ttk.Combobox(newwindow,state="readonly",values=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    combobox14.place(width=80,height=20,x=130,y=140)
    combobox14.current(int(alpha_plot_list_number)) 

    color_bar_list_0=["red","blue","green","yellow","gray"]
    if color_bar in color_bar_list_0:
        color_bar_list = ["red","blue","green","yellow","gray"]
    else:
        color_bar_list = ["red","blue","green","yellow","gray",color_bar]
    color_bar_list_number=color_bar_list.index(color_bar)
    global combobox15
    combobox15=ttk.Combobox(newwindow,state="readonly",values=color_bar_list)
    combobox15.place(width=80,height=20,x=130,y=200)
    combobox15.current(int(color_bar_list_number))  

    def color_choose_bar():
        color_now = color_bar
        color = colorchooser.askcolor(color_now,title="その他の色",master=newwindow)
        color_list=list(color)
        colorcode=color_list[1]
        if not colorcode == "None":
            global combobox15
            combobox15=ttk.Combobox(newwindow,state="readonly",values=["red","blue","green","yellow","gray",colorcode])
            combobox15.place(width=80,height=20,x=130,y=200)
            combobox15.current(5)  
        newwindow.attributes("-topmost",True)

    btn13=tk.Button(newwindow,text="その他の色",command=color_choose_bar)
    btn13.place(width=80,height=25,x=225,y=200)          

    bln2 = tkinter.BooleanVar()
    if bar_style == "枠あり":
        bln2.set(True)
    elif bar_style == "枠なし":
        bln2.set(False)
    else:
        print("checkbox_error")   
    check2 = tk.Checkbutton(newwindow,variable=bln2,text="棒グラフに枠線を表示する",font=(label_font_jp,9))
    check2.place(x=20,y=50)

    textBox34=ttk.Entry(newwindow)
    textBox34.place(width=60,height=20,x=130,y=110)
    textBox34.insert(tkinter.END,size_plot)

    bln3 = tkinter.BooleanVar()
    if tri_plot_color == "ON":
        bln3.set(True)
    elif tri_plot_color == "OFF":
        bln3.set(False)
    else:
        print("checkbox_error")   
    check3 = tk.Checkbutton(newwindow,variable=bln3,text="三角図上に分布を表示する",font=(label_font_jp,9))
    check3.place(x=20,y=230)

    def OK():
        global color_fig
        global color_axes
        global graph_font
        global color_plot
        global alpha_plot
        global color_bar
        global bar_style
        global size_plot
        global tri_plot_color
        if bln1.get():
            color_fig="black"
            color_axes="white"
        else:
            color_fig="white"
            color_axes="black"            
        print("color_fig:"+color_fig)
        print("color_axes:"+color_axes)
        graph_font=combobox12.get()
        print("graph_font:"+graph_font)
        color_plot=combobox13.get()
        print("color_plot:"+color_plot)
        alpha_plot=float(combobox14.get())
        print("alpha_plot:"+str(alpha_plot))
        color_bar=combobox15.get()
        print("color_bar:"+color_bar)    
        if bln2.get():
            bar_style="枠あり"
        else:
            bar_style="枠なし"
        print("bar_style:"+bar_style)  
        size_plot=int(float(textBox34.get()))
        print("size_plot:"+str(size_plot))  
        if bln3.get():
            tri_plot_color = "ON"
        else:
            tri_plot_color = "OFF"
        print("tri_plot_color:"+tri_plot_color)                   
        newwindow.destroy()

    def Cancel():
        global color_fig
        global color_axes
        global graph_font
        global color_plot
        global alpha_plot
        global color_bar
        global bar_style
        global size_plot
        global tri_plot_color
        color_fig=color_fig
        print("color_fig:"+color_fig)
        color_axes=color_axes
        print("color_axes:"+color_axes)
        graph_font=graph_font
        print("graph_font:"+graph_font)
        color_plot=color_plot
        print("color_plot:"+color_plot)
        alpha_plot=float(alpha_plot)
        print("alpha_plot:"+str(alpha_plot))
        color_bar=color_bar
        print("color_bar:"+color_bar)
        bar_style=bar_style
        print("bar_style:"+bar_style)
        size_plot=float(size_plot)
        print("size_plot:"+str(size_plot))
        tri_plot_color = tri_plot_color
        print("tri_plot_color:"+tri_plot_color)                   
        newwindow.destroy()

    btn9=tk.Button(newwindow,text="決定",command=OK)
    btn9.place(width=80,height=25,x=20,y=270)
    btn10=tk.Button(newwindow,text="キャンセル",command=Cancel)
    btn10.place(width=80,height=25,x=110,y=270)

def distribution3d():
    newwindow5 = tk.Toplevel(root)
    newwindow5.geometry("600x300")
    newwindow5.title(u"三次元組成分布設定")

    frame10 = tk.Canvas(newwindow5,width=80,height=30,bg="red")
    frame10.place(x=30,y=30)
    frame11 = tk.Canvas(newwindow5,width=80,height=30,bg="orange")
    frame11.place(x=110,y=30)
    frame12 = tk.Canvas(newwindow5,width=80,height=30,bg="lawngreen")
    frame12.place(x=190,y=30)
    frame13 = tk.Canvas(newwindow5,width=80,height=30,bg="green")
    frame13.place(x=270,y=30)
    frame14 = tk.Canvas(newwindow5,width=80,height=30,bg="blue")
    frame14.place(x=350,y=30)
    frame15 = tk.Canvas(newwindow5,width=80,height=30,bg="lightskyblue")
    frame15.place(x=430,y=30)

    textBox59 = ttk.Entry(newwindow5)
    textBox59.place(width=50,height=20,x=85,y=70)
    textBox59.insert(tkinter.END,hi_bord_1)
    textBox60 = ttk.Entry(newwindow5)
    textBox60.place(width=50,height=20,x=165,y=70)
    textBox60.insert(tkinter.END,hi_bord_2)
    textBox61 = ttk.Entry(newwindow5)
    textBox61.place(width=50,height=20,x=245,y=70)
    textBox61.insert(tkinter.END,hi_bord_3)
    textBox62 = ttk.Entry(newwindow5)
    textBox62.place(width=50,height=20,x=325,y=70)
    textBox62.insert(tkinter.END,hi_bord_4)
    textBox63 = ttk.Entry(newwindow5)
    textBox63.place(width=50,height=20,x=405,y=70)
    textBox63.insert(tkinter.END,hi_bord_5)

    def OK5():
        global hi_bord_1
        global hi_bord_2
        global hi_bord_3
        global hi_bord_4
        global hi_bord_5

        if int(textBox59.get()) >= int(textBox60.get()) >= int(textBox61.get()) >= int(textBox62.get()) >= int(textBox63.get()):
            hi_bord_1 = int(textBox59.get())
            hi_bord_2 = int(textBox60.get())
            hi_bord_3 = int(textBox61.get())
            hi_bord_4 = int(textBox62.get())
            hi_bord_5 = int(textBox63.get())

            newwindow5.destroy()
        else:
            messagebox.showerror("エラー","大小関係が不適切です")

    def Cancel5():
        global hi_bord_1
        global hi_bord_2
        global hi_bord_3
        global hi_bord_4
        global hi_bord_5

        hi_bord_1 = hi_bord_1
        hi_bord_2 = hi_bord_2
        hi_bord_3 = hi_bord_3
        hi_bord_4 = hi_bord_4
        hi_bord_5 = hi_bord_5

        newwindow5.destroy()

    btn20 = tk.Button(newwindow5,text="決定",command=OK5)
    btn20.place(width=80,height=25,x=400,y=120)
    btn21 = tk.Button(newwindow5,text="キャンセル",command=Cancel5)
    btn21.place(width=80,height=25,x=490,y=120)

def count_size():
    newwindow7 = tk.Toplevel(root)
    newwindow7.geometry("250x200")
    newwindow7.title(u"カウント数 → 粒径")

    label68 = tk.Label(newwindow7,text="Size Factor = ",font=(label_font_jp,9),anchor=tk.W)
    label68.place(width=100,height=20,x=20,y=30)

    textBox64 = ttk.Entry(newwindow7)
    textBox64.place(width=50,height=20,x=120,y=30)
    textBox64.insert(tkinter.END,size_factor)    

    bln3 = tkinter.BooleanVar()
    if tate_iso == "A":
        bln3.set(True)
    elif tate_iso == "all":
        bln3.set(False)
    else:
        print("checkbox_error")   
    check3 = tk.Checkbutton(newwindow7,variable=bln3,text="{}の大きさのみにする".format(combobox5.get()),font=(label_font_jp,9))
    check3.place(x=20,y=80)

    def OK7():
        global size_hosei
        global size_factor
        global tate_iso
        if textBox64.get() == "":
            messagebox.showerror("エラー","Size Factorを入力してください")
        elif float(textBox64.get()) <= 0:
            messagebox.showerror("エラー","Size Factorは正の値を入力してください")
        else:
            size_hosei = "sizefactor"
            size_factor = float(textBox64.get())
            if bln3.get():
                tate_iso = "A"
            else:
                tate_iso = "all"      
            newwindow7.destroy()

    def Cancel7():
        global size_hosei
        size_hosei = ""
        newwindow7.destroy()

    btn22=tk.Button(newwindow7,text="決定",command=OK7)
    btn22.place(width=80,height=25,x=20,y=140)
    btn23=tk.Button(newwindow7,text="キャンセル",command=Cancel7)
    btn23.place(width=80,height=25,x=110,y=140)

def end():
    ret = messagebox.askyesno("確認","ウィンドウを閉じますか？")
    if ret == True:
        root.quit()
        sys.exit()

def isospe():
    data_skip_column = int(data_large.columns.get_loc("skip"))
    element_number = int(data_skip_column)-2
    data_large_legend = data_large.columns[1:data_skip_column-1]
    if element_number <= 19:
        newwindow6 = tk.Toplevel(root)
        newwindow6_hight = element_number*30+90
        newwindow6.geometry("300x"+str(newwindow6_hight)+"")
        bln_nw6 = {}
        for i in range(0,element_number):
            bln_nw6[i] = tkinter.BooleanVar()
            bln_nw6[i].set(False)
            check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
            check_nw6_i.place(x=20,y=i*30+20)
    elif element_number <= 39:
        newwindow6 = tk.Toplevel(root)
        newwindow6_hight = 20*30+90
        newwindow6.geometry("600x"+str(newwindow6_hight)+"")
        bln_nw6 = {}
        for i in range(0,element_number):
            if i <= 19:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=20,y=i*30+20)
            elif i <= 39:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=320,y=(i-20)*30+20)
    elif element_number <= 59:
        newwindow6 = tk.Toplevel(root)
        newwindow6_hight = 20*30+90
        newwindow6.geometry("900x"+str(newwindow6_hight)+"")
        bln_nw6 = {}
        for i in range(0,element_number):
            if i <= 19:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=20,y=i*30+20)
            elif i <= 39:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=320,y=(i-20)*30+20)
            elif i <= 59:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=620,y=(i-40)*30+20)
    elif element_number <= 79:
        newwindow6 = tk.Toplevel(root)
        newwindow6_hight = 20*30+90
        newwindow6.geometry("1200x"+str(newwindow6_hight)+"")
        bln_nw6 = {}
        for i in range(0,element_number):
            if i <= 19:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=20,y=i*30+20)
            elif i <= 39:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=320,y=(i-20)*30+20)
            elif i <= 59:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=620,y=(i-40)*30+20)
            elif i <= 79:
                bln_nw6[i] = tkinter.BooleanVar()
                bln_nw6[i].set(False)
                check_nw6_i = tk.Checkbutton(newwindow6,variable=bln_nw6[i],text=data_large_legend[i],font=(label_font_jp,9))
                check_nw6_i.place(x=920,y=(i-60)*30+20)
    else:
        messagebox.showerror("エラー","測定同位体数が多すぎます")
    
    def OK6():
        global iso_spe
        global data_large_legend_spe
        global data_large_legend_spe_index
        global element_number_spe
        iso_spe = "specified"
        data_large_legend_spe = []
        data_large_legend_spe_index = []
        for i in range(0,element_number):
            if bln_nw6[i].get():
                data_large_legend_spe.append(data_large_legend[i])
                data_large_legend_spe_index.append(i)
        element_number_spe = len(data_large_legend_spe)
        newwindow6.destroy()

    def Cancel6():
        global iso_spe
        iso_spe = ""
        newwindow6.destroy()

    btn18 = tk.Button(newwindow6,text="決定",command=OK6)
    btn18.place(width=80,height=25,x=20,y=newwindow6_hight-50)
    btn19=tk.Button(newwindow6,text="キャンセル",command=Cancel6)
    btn19.place(width=80,height=25,x=110,y=newwindow6_hight-50)

def sansyou():
    newwindow2 = tk.Toplevel(root)
    newwindow2.geometry("500x540")
    newwindow2.title(u"参照値")

    label38=tk.Label(newwindow2,text="参照値１",font=(label_font_jp,9),anchor=tk.W)
    label38.place(width=60,height=20,x=20,y=20)
    if combobox5.get() == "":
        label39=tk.Label(newwindow2,text="元素A",font=(label_font_jp,9),anchor=tk.W)
    else:
        label39=tk.Label(newwindow2,text=combobox5.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label39.place(width=60,height=20,x=80,y=20)
    if combobox6.get() =="":
        label40=tk.Label(newwindow2,text="元素B",font=(label_font_jp,9),anchor=tk.W)
    else:
        label40=tk.Label(newwindow2,text=combobox6.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label40.place(width=60,height=20,x=80,y=50)
    if combobox7.get() =="":
        label41=tk.Label(newwindow2,text="元素C",font=(label_font_jp,9),anchor=tk.W)
    else:
        label41=tk.Label(newwindow2,text=combobox7.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label41.place(width=60,height=20,x=80,y=80)

    label42=tk.Label(newwindow2,text="参照値２",font=(label_font_jp,9),anchor=tk.W)
    label42.place(width=60,height=20,x=20,y=140)
    if combobox5.get() == "":
        label43=tk.Label(newwindow2,text="元素A",font=(label_font_jp,9),anchor=tk.W)
    else:
        label43=tk.Label(newwindow2,text=combobox5.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label43.place(width=60,height=20,x=80,y=140)
    if combobox6.get() =="":
        label44=tk.Label(newwindow2,text="元素B",font=(label_font_jp,9),anchor=tk.W)
    else:
        label44=tk.Label(newwindow2,text=combobox6.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label44.place(width=60,height=20,x=80,y=170)
    if combobox7.get() =="":
        label45=tk.Label(newwindow2,text="元素C",font=(label_font_jp,9),anchor=tk.W)
    else:
        label45=tk.Label(newwindow2,text=combobox7.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label45.place(width=60,height=20,x=80,y=200)

    label46=tk.Label(newwindow2,text="参照値３",font=(label_font_jp,9),anchor=tk.W)
    label46.place(width=60,height=20,x=20,y=260)
    if combobox5.get() == "":
        label47=tk.Label(newwindow2,text="元素A",font=(label_font_jp,9),anchor=tk.W)
    else:
        label47=tk.Label(newwindow2,text=combobox5.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label47.place(width=60,height=20,x=80,y=260)
    if combobox6.get() =="":
        label48=tk.Label(newwindow2,text="元素B",font=(label_font_jp,9),anchor=tk.W)
    else:
        label48=tk.Label(newwindow2,text=combobox6.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label48.place(width=60,height=20,x=80,y=290)
    if combobox7.get() =="":
        label49=tk.Label(newwindow2,text="元素C",font=(label_font_jp,9),anchor=tk.W)
    else:
        label49=tk.Label(newwindow2,text=combobox7.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label49.place(width=60,height=20,x=80,y=320)

    label46_4=tk.Label(newwindow2,text="参照値４",font=(label_font_jp,9),anchor=tk.W)
    label46_4.place(width=60,height=20,x=20,y=380)
    if combobox5.get() == "":
        label47_4=tk.Label(newwindow2,text="元素A",font=(label_font_jp,9),anchor=tk.W)
    else:
        label47_4=tk.Label(newwindow2,text=combobox5.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label47_4.place(width=60,height=20,x=80,y=380)
    if combobox6.get() =="":
        label48_4=tk.Label(newwindow2,text="元素B",font=(label_font_jp,9),anchor=tk.W)
    else:
        label48_4=tk.Label(newwindow2,text=combobox6.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label48_4.place(width=60,height=20,x=80,y=410)
    if combobox7.get() =="":
        label49_4=tk.Label(newwindow2,text="元素C",font=(label_font_jp,9),anchor=tk.W)
    else:
        label49_4=tk.Label(newwindow2,text=combobox7.get(),font=(label_font_jp,9),anchor=tk.CENTER)
    label49_4.place(width=60,height=20,x=80,y=440)

    textBox35=ttk.Entry(newwindow2)
    textBox35.place(width=50,height=20,x=140,y=20)
    textBox35.insert(tkinter.END,sansyou_1_a)
    label50=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label50.place(width=20,height=20,x=190,y=20)
    textBox36=ttk.Entry(newwindow2)
    textBox36.place(width=50,height=20,x=210,y=20)
    textBox36.insert(tkinter.END,sansyou_1_sa)
    textBox37=ttk.Entry(newwindow2)
    textBox37.place(width=50,height=20,x=140,y=50)
    textBox37.insert(tkinter.END,sansyou_1_b)
    label51=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label51.place(width=20,height=20,x=190,y=50)
    textBox38=ttk.Entry(newwindow2)
    textBox38.place(width=50,height=20,x=210,y=50)
    textBox38.insert(tkinter.END,sansyou_1_sb)
    textBox39=ttk.Entry(newwindow2)
    textBox39.place(width=50,height=20,x=140,y=80)
    textBox39.insert(tkinter.END,sansyou_1_c)
    label52=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label52.place(width=20,height=20,x=190,y=80)
    textBox40=ttk.Entry(newwindow2)
    textBox40.place(width=50,height=20,x=210,y=80)
    textBox40.insert(tkinter.END,sansyou_1_sc)

    textBox41=ttk.Entry(newwindow2)
    textBox41.place(width=50,height=20,x=140,y=140)
    textBox41.insert(tkinter.END,sansyou_2_a)
    label53=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label53.place(width=20,height=20,x=190,y=140)
    textBox42=ttk.Entry(newwindow2)
    textBox42.place(width=50,height=20,x=210,y=140)
    textBox42.insert(tkinter.END,sansyou_2_sa)
    textBox43=ttk.Entry(newwindow2)
    textBox43.place(width=50,height=20,x=140,y=170)
    textBox43.insert(tkinter.END,sansyou_2_b)
    label54=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label54.place(width=20,height=20,x=190,y=170)
    textBox44=ttk.Entry(newwindow2)
    textBox44.place(width=50,height=20,x=210,y=170)
    textBox44.insert(tkinter.END,sansyou_2_sb)
    textBox45=ttk.Entry(newwindow2)
    textBox45.place(width=50,height=20,x=140,y=200)
    textBox45.insert(tkinter.END,sansyou_2_c)
    label55=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label55.place(width=20,height=20,x=190,y=200)
    textBox46=ttk.Entry(newwindow2)
    textBox46.place(width=50,height=20,x=210,y=200)
    textBox46.insert(tkinter.END,sansyou_2_sc)

    textBox47=ttk.Entry(newwindow2)
    textBox47.place(width=50,height=20,x=140,y=260)
    textBox47.insert(tkinter.END,sansyou_3_a)
    label56=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label56.place(width=20,height=20,x=190,y=260)
    textBox48=ttk.Entry(newwindow2)
    textBox48.place(width=50,height=20,x=210,y=260)
    textBox48.insert(tkinter.END,sansyou_3_sa)
    textBox49=ttk.Entry(newwindow2)
    textBox49.place(width=50,height=20,x=140,y=290)
    textBox49.insert(tkinter.END,sansyou_3_b)
    label57=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label57.place(width=20,height=20,x=190,y=290)
    textBox50=ttk.Entry(newwindow2)
    textBox50.place(width=50,height=20,x=210,y=290)
    textBox50.insert(tkinter.END,sansyou_3_sb)
    textBox51=ttk.Entry(newwindow2)
    textBox51.place(width=50,height=20,x=140,y=320)
    textBox51.insert(tkinter.END,sansyou_3_c)
    label58=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label58.place(width=20,height=20,x=190,y=320)
    textBox52=ttk.Entry(newwindow2)
    textBox52.place(width=50,height=20,x=210,y=320)
    textBox52.insert(tkinter.END,sansyou_3_sc)

    textBox47_4=ttk.Entry(newwindow2)
    textBox47_4.place(width=50,height=20,x=140,y=380)
    textBox47_4.insert(tkinter.END,sansyou_4_a)
    label56_4=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label56_4.place(width=20,height=20,x=190,y=380)
    textBox48_4=ttk.Entry(newwindow2)
    textBox48_4.place(width=50,height=20,x=210,y=380)
    textBox48_4.insert(tkinter.END,sansyou_4_sa)
    textBox49_4=ttk.Entry(newwindow2)
    textBox49_4.place(width=50,height=20,x=140,y=410)
    textBox49_4.insert(tkinter.END,sansyou_4_b)
    label57_4=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label57_4.place(width=20,height=20,x=190,y=410)
    textBox50_4=ttk.Entry(newwindow2)
    textBox50_4.place(width=50,height=20,x=210,y=410)
    textBox50_4.insert(tkinter.END,sansyou_4_sb)
    textBox51_4=ttk.Entry(newwindow2)
    textBox51_4.place(width=50,height=20,x=140,y=440)
    textBox51_4.insert(tkinter.END,sansyou_4_c)
    label58_4=tk.Label(newwindow2,text="±",font=(label_font_jp,9),anchor=tk.CENTER)
    label58_4.place(width=20,height=20,x=190,y=440)
    textBox52_4=ttk.Entry(newwindow2)
    textBox52_4.place(width=50,height=20,x=210,y=440)
    textBox52_4.insert(tkinter.END,sansyou_4_sc)

    def OK2():
        global sansyou_1_a
        global sansyou_1_b
        global sansyou_1_c
        global sansyou_1_sa
        global sansyou_1_sb
        global sansyou_1_sc
        global sansyou_2_a
        global sansyou_2_b
        global sansyou_2_c
        global sansyou_2_sa
        global sansyou_2_sb
        global sansyou_2_sc
        global sansyou_3_a
        global sansyou_3_b
        global sansyou_3_c
        global sansyou_3_sa
        global sansyou_3_sb
        global sansyou_3_sc
        global sansyou_4_a
        global sansyou_4_b
        global sansyou_4_c
        global sansyou_4_sa
        global sansyou_4_sb
        global sansyou_4_sc

        sansyou_1_a = textBox35.get()
        sansyou_1_sa = textBox36.get()
        sansyou_1_b = textBox37.get()
        sansyou_1_sb = textBox38.get()
        sansyou_1_c = textBox39.get()
        sansyou_1_sc = textBox40.get()

        sansyou_2_a = textBox41.get()
        sansyou_2_sa = textBox42.get()
        sansyou_2_b = textBox43.get()
        sansyou_2_sb = textBox44.get()
        sansyou_2_c = textBox45.get()
        sansyou_2_sc = textBox46.get()

        sansyou_3_a = textBox47.get()
        sansyou_3_sa = textBox48.get()
        sansyou_3_b = textBox49.get()
        sansyou_3_sb = textBox50.get()
        sansyou_3_c = textBox51.get()
        sansyou_3_sc = textBox52.get()

        sansyou_4_a = textBox47_4.get()
        sansyou_4_sa = textBox48_4.get()
        sansyou_4_b = textBox49_4.get()
        sansyou_4_sb = textBox50_4.get()
        sansyou_4_c = textBox51_4.get()
        sansyou_4_sc = textBox52_4.get()

        newwindow2.destroy()

    def Cancel2():
        global sansyou_1_a
        global sansyou_1_b
        global sansyou_1_c
        global sansyou_1_sa
        global sansyou_1_sb
        global sansyou_1_sc
        global sansyou_2_a
        global sansyou_2_b
        global sansyou_2_c
        global sansyou_2_sa
        global sansyou_2_sb
        global sansyou_2_sc
        global sansyou_3_a
        global sansyou_3_b
        global sansyou_3_c
        global sansyou_3_sa
        global sansyou_3_sb
        global sansyou_3_sc
        global sansyou_4_a
        global sansyou_4_b
        global sansyou_4_c
        global sansyou_4_sa
        global sansyou_4_sb
        global sansyou_4_sc


        sansyou_1_a = sansyou_1_a
        sansyou_1_b = sansyou_1_b
        sansyou_1_c = sansyou_1_c
        sansyou_1_sa = sansyou_1_sa
        sansyou_1_sb = sansyou_1_sb
        sansyou_1_sc = sansyou_1_sc

        sansyou_2_a = sansyou_2_a
        sansyou_2_b = sansyou_2_b
        sansyou_2_c = sansyou_2_c
        sansyou_2_sa = sansyou_2_sa
        sansyou_2_sb = sansyou_2_sb
        sansyou_2_sc = sansyou_2_sc

        sansyou_3_a = sansyou_3_a
        sansyou_3_b = sansyou_3_b
        sansyou_3_c = sansyou_3_c
        sansyou_3_sa = sansyou_3_sa
        sansyou_3_sb = sansyou_3_sb
        sansyou_3_sc = sansyou_3_sc

        sansyou_4_a = sansyou_4_a
        sansyou_4_b = sansyou_4_b
        sansyou_4_c = sansyou_4_c
        sansyou_4_sa = sansyou_4_sa
        sansyou_4_sb = sansyou_4_sb
        sansyou_4_sc = sansyou_4_sc

        newwindow2.destroy()
    
    btn13=tk.Button(newwindow2,text="決定",command=OK2)
    btn13.place(width=80,height=25,x=20,y=490)
    btn14=tk.Button(newwindow2,text="キャンセル",command=Cancel2)
    btn14.place(width=80,height=25,x=110,y=490)

def choutenplus():
    if textBox2.get() == "":
        messagebox.showerror("エラー","ファイルを選択してください")
    else:
        data_csv_p = pd.read_csv(textBox2.get())
        isotope_df_0_p = data_csv_p.columns
        isonum = (len(isotope_df_0_p)-4)/2
        isotope_df_p = isotope_df_0_p[1:int(isonum)+1]
        isotope_list_p = isotope_df_p.values.tolist()
        ele_list_p = []
        for i in range(0,len(isotope_list_p)):
            alpha_p_i = "".join([s for s in isotope_list_p[i] if s.isalpha()])
            ele_list_p.append(alpha_p_i)
        mass_list_p = []
        for i in range(0,len(isotope_list_p)):
            digit_p_i = "".join([s for s in isotope_list_p[i] if s.isdigit()])
            mass_list_p.append(digit_p_i)
        ele_list_min_p = list(set(ele_list_p))
        ele_list_min_p.sort()
        ele_list_min_ini_p = list(set(ele_list_p))
        ele_list_min_ini_p.sort()
        ele_list_min_bar_p = list(set(ele_list_p))
        ele_list_min_bar_p.sort()
        ele_list_min_bar_p.insert(0,"-")
        ele_list_min_zero_p = list(set(ele_list_p))
        ele_list_min_zero_p.sort()
        ele_list_min_zero_p.append("-")

        if len(ele_list_min_p) > 15:
            messagebox.showerror("エラー","元素数が多すぎます")
        else:
            while len(ele_list_min_p) < 15:
                ele_list_min_p.append('-')
        
        ele_list_use_p = []

        if ele_list_min_p[0] == "-":
            messagebox.showerror('エラー',"最低３元素は測定してください")
        else:
            ele_list_use_p.append(ele_list_min_p[0])
            iso_list_p1 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[0]]
            mass_list_p1 = []
            for i in range(0,len(iso_list_p1)):
                iso_index_p_i = iso_list_p1[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p1.append(mass_p_i)
        if ele_list_min_p[1] == "-":
            messagebox.showerror('エラー',"最低３元素は測定してください")
        else:
            ele_list_use_p.append(ele_list_min_p[1])
            iso_list_p2 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[1]]
            mass_list_p2 = []
            for i in range(0,len(iso_list_p2)):
                iso_index_p_i = iso_list_p2[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p2.append(mass_p_i)
        if ele_list_min_p[2] == "-":
            messagebox.showerror('エラー',"最低３元素は測定してください")
        else:
            ele_list_use_p.append(ele_list_min_p[2])
            iso_list_p3 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[2]]
            mass_list_p3 = []
            for i in range(0,len(iso_list_p3)):
                iso_index_p_i = iso_list_p3[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p3.append(mass_p_i)
        if ele_list_min_p[3] == "-":
            mass_list_p4 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[3])
            iso_list_p4 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[3]]
            mass_list_p4 = []
            for i in range(0,len(iso_list_p4)):
                iso_index_p_i = iso_list_p4[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p4.append(mass_p_i)
        if ele_list_min_p[4] == "-":
            mass_list_p5 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[4])
            iso_list_p5 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[4]]
            mass_list_p5 = []
            for i in range(0,len(iso_list_p5)):
                iso_index_p_i = iso_list_p5[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p5.append(mass_p_i)
        if ele_list_min_p[5] == "-":
            mass_list_p6 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[5])
            iso_list_p6 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[5]]
            mass_list_p6 = []
            for i in range(0,len(iso_list_p6)):
                iso_index_p_i = iso_list_p6[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p6.append(mass_p_i)
        if ele_list_min_p[6] == "-":
            mass_list_p7 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[6])
            iso_list_p7 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[6]]
            mass_list_p7 = []
            for i in range(0,len(iso_list_p7)):
                iso_index_p_i = iso_list_p7[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p7.append(mass_p_i)
        if ele_list_min_p[7] == "-":
            mass_list_p8 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[7])
            iso_list_p8 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[7]]
            mass_list_p8 = []
            for i in range(0,len(iso_list_p8)):
                iso_index_p_i = iso_list_p8[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p8.append(mass_p_i)
        if ele_list_min_p[8] == "-":
            mass_list_p9 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[8])
            iso_list_p9 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[8]]
            mass_list_p9 = []
            for i in range(0,len(iso_list_p9)):
                iso_index_p_i = iso_list_p9[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p9.append(mass_p_i)
        if ele_list_min_p[9] == "-":
            mass_list_p10 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[9])
            iso_list_p10 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[9]]
            mass_list_p10 = []
            for i in range(0,len(iso_list_p10)):
                iso_index_p_i = iso_list_p10[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p10.append(mass_p_i)
        if ele_list_min_p[10] == "-":
            mass_list_p11 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[10])
            iso_list_p11 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[10]]
            mass_list_p11 = []
            for i in range(0,len(iso_list_p11)):
                iso_index_p_i = iso_list_p11[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p11.append(mass_p_i)
        if ele_list_min_p[11] == "-":
            mass_list_p12 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[11])
            iso_list_p12 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[11]]
            mass_list_p12 = []
            for i in range(0,len(iso_list_p12)):
                iso_index_p_i = iso_list_p12[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p12.append(mass_p_i)
        if ele_list_min_p[12] == "-":
            mass_list_p13 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[12])
            iso_list_p13 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[12]]
            mass_list_p13 = []
            for i in range(0,len(iso_list_p13)):
                iso_index_p_i = iso_list_p13[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p13.append(mass_p_i)
        if ele_list_min_p[13] == "-":
            mass_list_p14 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[13])
            iso_list_p14 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[13]]
            mass_list_p14 = []
            for i in range(0,len(iso_list_p14)):
                iso_index_p_i = iso_list_p14[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p14.append(mass_p_i)
        if ele_list_min_p[14] == "-":
            mass_list_p15 = ["-"]
        else:
            ele_list_use_p.append(ele_list_min_p[14])
            iso_list_p15 = [k for k, x in enumerate(ele_list_p) if x == ele_list_min_p[14]]
            mass_list_p15 = []
            for i in range(0,len(iso_list_p15)):
                iso_index_p_i = iso_list_p15[i]
                mass_p_i = mass_list_p[int(iso_index_p_i)]
                mass_list_p15.append(mass_p_i)
        
        newwindow_p = tk.Toplevel(root)
        newwindow_p.attributes("-topmost",True)
        newwindow_p.geometry("1290x280")
        newwindow_p.title(u"元素読み込み")

        label_p01 = tk.Label(newwindow_p,text="元素",font=(label_font_jp,10),anchor=tk.CENTER)
        label_p01.place(width=50,height=20,x=20,y=20)
        label_p02 = tk.Label(newwindow_p,text="質量数",font=(label_font_jp,10),anchor=tk.CENTER)
        label_p02.place(width=50,height=20,x=20,y=60)
        label_p03 = tk.Label(newwindow_p,text="補正",font=(label_font_jp,10),anchor=tk.CENTER)
        label_p03.place(width=50,height=20,x=20,y=100)
        label_p03 = tk.Label(newwindow_p,text="サイズ",font=(label_font_jp,10),anchor=tk.CENTER)
        label_p03.place(width=50,height=20,x=20,y=140)
        label_p03 = tk.Label(newwindow_p,text="基準",font=(label_font_jp,10),anchor=tk.CENTER)
        label_p03.place(width=50,height=20,x=20,y=180)

        label_p1 = tk.Label(newwindow_p,text=ele_list_min_p[0],font=(label_font_jp,10))
        label_p1.place(width=50,height=20,x=90,y=20)
        label_p2 = tk.Label(newwindow_p,text=ele_list_min_p[1],font=(label_font_jp,10))
        label_p2.place(width=50,height=20,x=170,y=20)
        label_p3 = tk.Label(newwindow_p,text=ele_list_min_p[2],font=(label_font_jp,10))
        label_p3.place(width=50,height=20,x=250,y=20)
        label_p4 = tk.Label(newwindow_p,text=ele_list_min_p[3],font=(label_font_jp,10))
        label_p4.place(width=50,height=20,x=330,y=20)
        label_p5 = tk.Label(newwindow_p,text=ele_list_min_p[4],font=(label_font_jp,10))
        label_p5.place(width=50,height=20,x=410,y=20)
        label_p6 = tk.Label(newwindow_p,text=ele_list_min_p[5],font=(label_font_jp,10))
        label_p6.place(width=50,height=20,x=490,y=20)
        label_p7 = tk.Label(newwindow_p,text=ele_list_min_p[6],font=(label_font_jp,10))
        label_p7.place(width=50,height=20,x=570,y=20)
        label_p8 = tk.Label(newwindow_p,text=ele_list_min_p[7],font=(label_font_jp,10))
        label_p8.place(width=50,height=20,x=650,y=20)
        label_p9 = tk.Label(newwindow_p,text=ele_list_min_p[8],font=(label_font_jp,10))
        label_p9.place(width=50,height=20,x=730,y=20)
        label_p10 = tk.Label(newwindow_p,text=ele_list_min_p[9],font=(label_font_jp,10))
        label_p10.place(width=50,height=20,x=810,y=20)
        label_p011 = tk.Label(newwindow_p,text=ele_list_min_p[10],font=(label_font_jp,10))
        label_p011.place(width=50,height=20,x=890,y=20)
        label_p012 = tk.Label(newwindow_p,text=ele_list_min_p[11],font=(label_font_jp,10))
        label_p012.place(width=50,height=20,x=970,y=20)
        label_p013 = tk.Label(newwindow_p,text=ele_list_min_p[12],font=(label_font_jp,10))
        label_p013.place(width=50,height=20,x=1050,y=20)
        label_p014 = tk.Label(newwindow_p,text=ele_list_min_p[13],font=(label_font_jp,10))
        label_p014.place(width=50,height=20,x=1130,y=20)
        label_p015 = tk.Label(newwindow_p,text=ele_list_min_p[14],font=(label_font_jp,10))
        label_p015.place(width=50,height=20,x=1210,y=20)

        combobox_p1 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p1)
        combobox_p1.place(width=50,height=20,x=90,y=60)
        combobox_p1.current(0)
        combobox_p2 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p2)
        combobox_p2.place(width=50,height=20,x=170,y=60)
        combobox_p2.current(0)
        combobox_p3 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p3)
        combobox_p3.place(width=50,height=20,x=250,y=60)
        combobox_p3.current(0)
        combobox_p4 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p4)
        combobox_p4.place(width=50,height=20,x=330,y=60)
        combobox_p4.current(0)
        combobox_p5 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p5)
        combobox_p5.place(width=50,height=20,x=410,y=60)
        combobox_p5.current(0)
        combobox_p6 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p6)
        combobox_p6.place(width=50,height=20,x=490,y=60)
        combobox_p6.current(0)
        combobox_p7 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p7)
        combobox_p7.place(width=50,height=20,x=570,y=60)
        combobox_p7.current(0)
        combobox_p8 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p8)
        combobox_p8.place(width=50,height=20,x=650,y=60)
        combobox_p8.current(0)
        combobox_p9 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p9)
        combobox_p9.place(width=50,height=20,x=730,y=60)
        combobox_p9.current(0)
        combobox_p10 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p10)
        combobox_p10.place(width=50,height=20,x=810,y=60)
        combobox_p10.current(0)
        combobox_p011 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p11)
        combobox_p011.place(width=50,height=20,x=890,y=60)
        combobox_p011.current(0)
        combobox_p012 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p12)
        combobox_p012.place(width=50,height=20,x=970,y=60)
        combobox_p012.current(0)
        combobox_p013 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p13)
        combobox_p013.place(width=50,height=20,x=1050,y=60)
        combobox_p013.current(0)
        combobox_p014 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p14)
        combobox_p014.place(width=50,height=20,x=1130,y=60)
        combobox_p014.current(0)
        combobox_p015 = ttk.Combobox(newwindow_p,state='readonly',values=mass_list_p15)
        combobox_p015.place(width=50,height=20,x=1210,y=60)
        combobox_p015.current(0)

        textBox_p1 = ttk.Entry(newwindow_p)
        textBox_p1.place(width=50,height=20,x=90,y=100)
        textBox_p1.insert(tkinter.END,cor_p1)
        textBox_p2 = ttk.Entry(newwindow_p)
        textBox_p2.place(width=50,height=20,x=170,y=100)
        textBox_p2.insert(tkinter.END,cor_p2)
        textBox_p3 = ttk.Entry(newwindow_p)
        textBox_p3.place(width=50,height=20,x=250,y=100)
        textBox_p3.insert(tkinter.END,cor_p3)
        textBox_p4 = ttk.Entry(newwindow_p)
        textBox_p4.place(width=50,height=20,x=330,y=100)
        textBox_p4.insert(tkinter.END,cor_p4)
        textBox_p5 = ttk.Entry(newwindow_p)
        textBox_p5.place(width=50,height=20,x=410,y=100)
        textBox_p5.insert(tkinter.END,cor_p5)
        textBox_p6 = ttk.Entry(newwindow_p)
        textBox_p6.place(width=50,height=20,x=490,y=100)
        textBox_p6.insert(tkinter.END,cor_p6)
        textBox_p7 = ttk.Entry(newwindow_p)
        textBox_p7.place(width=50,height=20,x=570,y=100)
        textBox_p7.insert(tkinter.END,cor_p7)
        textBox_p8 = ttk.Entry(newwindow_p)
        textBox_p8.place(width=50,height=20,x=650,y=100)
        textBox_p8.insert(tkinter.END,cor_p8)
        textBox_p9 = ttk.Entry(newwindow_p)
        textBox_p9.place(width=50,height=20,x=730,y=100)
        textBox_p9.insert(tkinter.END,cor_p9)
        textBox_p10 = ttk.Entry(newwindow_p)
        textBox_p10.place(width=50,height=20,x=810,y=100)
        textBox_p10.insert(tkinter.END,cor_p10)
        textBox_p011 = ttk.Entry(newwindow_p)
        textBox_p011.place(width=50,height=20,x=890,y=100)
        textBox_p011.insert(tkinter.END,cor_p11)
        textBox_p012 = ttk.Entry(newwindow_p)
        textBox_p012.place(width=50,height=20,x=970,y=100)
        textBox_p012.insert(tkinter.END,cor_p12)
        textBox_p013 = ttk.Entry(newwindow_p)
        textBox_p013.place(width=50,height=20,x=1050,y=100)
        textBox_p013.insert(tkinter.END,cor_p13)
        textBox_p014 = ttk.Entry(newwindow_p)
        textBox_p014.place(width=50,height=20,x=1130,y=100)
        textBox_p014.insert(tkinter.END,cor_p14)
        textBox_p015 = ttk.Entry(newwindow_p)
        textBox_p015.place(width=50,height=20,x=1210,y=100)
        textBox_p015.insert(tkinter.END,cor_p15)
        
        global bln_p1
        global bln_p2
        global bln_p3
        global bln_p4
        global bln_p5
        global bln_p6
        global bln_p7
        global bln_p8
        global bln_p9
        global bln_p10        
        global bln_p11  
        global bln_p12        
        global bln_p13        
        global bln_p14        
        global bln_p15        
        bln_p1 = tkinter.BooleanVar()
        if ele_list_min_p[0] == "-":
            bln_p1.set(False)      
        elif size_calc_p1 == "true":
            bln_p1.set(True)
        elif size_calc_p1 == "false":
            bln_p1.set(False)
        check_p1 = tk.Checkbutton(newwindow_p,variable=bln_p1)
        check_p1.place(x=90,y=140)
        bln_p2 = tkinter.BooleanVar()
        if ele_list_min_p[1] == "-":
            bln_p2.set(False)        
        elif size_calc_p2 == "true":
            bln_p2.set(True)
        elif size_calc_p2 == "false":
            bln_p2.set(False)
        else:
            print("checkbox_error")   
        check_p2 = tk.Checkbutton(newwindow_p,variable=bln_p2)
        check_p2.place(x=170,y=140)
        bln_p3 = tkinter.BooleanVar()
        if ele_list_min_p[2] == "-":
            bln_p3.set(False)        
        elif size_calc_p3 == "true":
            bln_p3.set(True)
        elif size_calc_p3 == "false":
            bln_p3.set(False)
        else:
            print("checkbox_error")   
        check_p3 = tk.Checkbutton(newwindow_p,variable=bln_p3)
        check_p3.place(x=250,y=140)
        bln_p4 = tkinter.BooleanVar()
        if ele_list_min_p[3] == "-":
            bln_p4.set(False)        
        elif size_calc_p4 == "true":
            bln_p4.set(True)
        elif size_calc_p4 == "false":
            bln_p4.set(False)
        else:
            print("checkbox_error")   
        check_p4 = tk.Checkbutton(newwindow_p,variable=bln_p4)
        check_p4.place(x=330,y=140)
        bln_p5 = tkinter.BooleanVar()
        if ele_list_min_p[4] == "-":
            bln_p5.set(False)        
        elif size_calc_p5 == "true":
            bln_p5.set(True)
        elif size_calc_p5 == "false":
            bln_p5.set(False)
        else:
            print("checkbox_error")   
        check_p5 = tk.Checkbutton(newwindow_p,variable=bln_p5)
        check_p5.place(x=410,y=140)
        bln_p6 = tkinter.BooleanVar()
        if ele_list_min_p[5] == "-":
            bln_p6.set(False)        
        elif size_calc_p6 == "true":
            bln_p6.set(True)
        elif size_calc_p6 == "false":
            bln_p6.set(False)
        else:
            print("checkbox_error")   
        check_p6 = tk.Checkbutton(newwindow_p,variable=bln_p6)
        check_p6.place(x=490,y=140)
        bln_p7 = tkinter.BooleanVar()
        if ele_list_min_p[6] == "-":
            bln_p7.set(False)        
        elif size_calc_p7 == "true":
            bln_p7.set(True)
        elif size_calc_p7 == "false":
            bln_p7.set(False)
        else:
            print("checkbox_error")   
        check_p7 = tk.Checkbutton(newwindow_p,variable=bln_p7)
        check_p7.place(x=570,y=140)
        bln_p8 = tkinter.BooleanVar()
        if ele_list_min_p[7] == "-":
            bln_p8.set(False)        
        elif size_calc_p8 == "true":
            bln_p8.set(True)
        elif size_calc_p8 == "false":
            bln_p8.set(False)
        else:
            print("checkbox_error")   
        check_p8 = tk.Checkbutton(newwindow_p,variable=bln_p8)
        check_p8.place(x=650,y=140)
        bln_p9 = tkinter.BooleanVar()
        if ele_list_min_p[8] == "-":
            bln_p9.set(False)        
        elif size_calc_p9 == "true":
            bln_p9.set(True)
        elif size_calc_p9 == "false":
            bln_p9.set(False)
        else:
            print("checkbox_error")   
        check_p9 = tk.Checkbutton(newwindow_p,variable=bln_p9)
        check_p9.place(x=730,y=140)
        bln_p10 = tkinter.BooleanVar()
        if ele_list_min_p[9] == "-":
            bln_p10.set(False)        
        elif size_calc_p10 == "true":
            bln_p10.set(True)
        elif size_calc_p10 == "false":
            bln_p10.set(False)
        else:
            print("checkbox_error")   
        check_p10 = tk.Checkbutton(newwindow_p,variable=bln_p10)
        check_p10.place(x=810,y=140)
        bln_p11 = tkinter.BooleanVar()
        if ele_list_min_p[10] == "-":
            bln_p11.set(False)        
        elif size_calc_p11 == "true":
            bln_p11.set(True)
        elif size_calc_p11 == "false":
            bln_p11.set(False)
        else:
            print("checkbox_error")   
        check_p11 = tk.Checkbutton(newwindow_p,variable=bln_p11)
        check_p11.place(x=890,y=140)
        bln_p12 = tkinter.BooleanVar()
        if ele_list_min_p[11] == "-":
            bln_p12.set(False)        
        elif size_calc_p12 == "true":
            bln_p12.set(True)
        elif size_calc_p12 == "false":
            bln_p12.set(False)
        else:
            print("checkbox_error")   
        check_p12 = tk.Checkbutton(newwindow_p,variable=bln_p12)
        check_p12.place(x=970,y=140)
        bln_p13 = tkinter.BooleanVar()
        if ele_list_min_p[12] == "-":
            bln_p13.set(False)        
        elif size_calc_p13 == "true":
            bln_p13.set(True)
        elif size_calc_p13 == "false":
            bln_p13.set(False)
        else:
            print("checkbox_error")   
        check_p13 = tk.Checkbutton(newwindow_p,variable=bln_p13)
        check_p13.place(x=1050,y=140)
        bln_p14 = tkinter.BooleanVar()
        if ele_list_min_p[13] == "-":
            bln_p14.set(False)        
        elif size_calc_p14 == "true":
            bln_p14.set(True)
        elif size_calc_p14 == "false":
            bln_p14.set(False)
        else:
            print("checkbox_error")   
        check_p14 = tk.Checkbutton(newwindow_p,variable=bln_p14)
        check_p14.place(x=1130,y=140)
        bln_p15 = tkinter.BooleanVar()
        if ele_list_min_p[14] == "-":
            bln_p15.set(False)        
        elif size_calc_p15 == "true":
            bln_p15.set(True)
        elif size_calc_p15 == "false":
            bln_p15.set(False)
        else:
            print("checkbox_error")   
        check_p15 = tk.Checkbutton(newwindow_p,variable=bln_p15)
        check_p15.place(x=1210,y=140)

        combobox_p23 = ttk.Combobox(newwindow_p,state='readonly',values=ele_list_use_p)
        combobox_p23.place(width=50,height=20,x=90,y=180)
        combobox_p23.current(0)

    def OKp():
        global ele_p1
        global ele_p2
        global ele_p3
        global ele_p4
        global ele_p5
        global ele_p6
        global ele_p7
        global ele_p8
        global ele_p9
        global ele_p10
        global ele_p11
        global ele_p12
        global ele_p13
        global ele_p14
        global ele_p15
        global mass_p1
        global mass_p2
        global mass_p3
        global mass_p4
        global mass_p5
        global mass_p6
        global mass_p7
        global mass_p8
        global mass_p9
        global mass_p10
        global mass_p11
        global mass_p12
        global mass_p13
        global mass_p14
        global mass_p15
        global cor_p1
        global cor_p2
        global cor_p3
        global cor_p4
        global cor_p5
        global cor_p6
        global cor_p7
        global cor_p8
        global cor_p9
        global cor_p10
        global cor_p11
        global cor_p12
        global cor_p13
        global cor_p14
        global cor_p15
        global mass_index_p_list

        ele_p1 = ele_list_min_p[0]
        ele_p2 = ele_list_min_p[1]
        ele_p3 = ele_list_min_p[2]
        ele_p4 = ele_list_min_p[3]
        ele_p5 = ele_list_min_p[4]
        ele_p6 = ele_list_min_p[5]
        ele_p7 = ele_list_min_p[6]
        ele_p8 = ele_list_min_p[7]
        ele_p9 = ele_list_min_p[8]
        ele_p10 = ele_list_min_p[9]
        ele_p11 = ele_list_min_p[10]
        ele_p12 = ele_list_min_p[11]
        ele_p13 = ele_list_min_p[12]
        ele_p14 = ele_list_min_p[13]
        ele_p15 = ele_list_min_p[14]
        ele_p_list = []
        ele_p_list.append(ele_p1)
        ele_p_list.append(ele_p2)
        ele_p_list.append(ele_p3)
        ele_p_list.append(ele_p4)
        ele_p_list.append(ele_p5)
        ele_p_list.append(ele_p6)
        ele_p_list.append(ele_p7)
        ele_p_list.append(ele_p8)
        ele_p_list.append(ele_p9)
        ele_p_list.append(ele_p10)
        ele_p_list.append(ele_p11)
        ele_p_list.append(ele_p12)
        ele_p_list.append(ele_p13)
        ele_p_list.append(ele_p14)
        ele_p_list.append(ele_p15)
        print(ele_p_list)

        mass_p1 = combobox_p1.get()
        mass_p2 = combobox_p2.get()
        mass_p3 = combobox_p3.get()
        mass_p4 = combobox_p4.get()
        mass_p5 = combobox_p5.get()
        mass_p6 = combobox_p6.get()
        mass_p7 = combobox_p7.get()
        mass_p8 = combobox_p8.get()
        mass_p9 = combobox_p9.get()
        mass_p10 = combobox_p10.get()
        mass_p11 = combobox_p011.get()
        mass_p12 = combobox_p012.get()
        mass_p13 = combobox_p013.get()
        mass_p14 = combobox_p014.get()
        mass_p15 = combobox_p015.get()
        mass_p_list =[]
        mass_p_list.append(mass_p1)
        mass_p_list.append(mass_p2)
        mass_p_list.append(mass_p3)
        mass_p_list.append(mass_p4)
        mass_p_list.append(mass_p5)
        mass_p_list.append(mass_p6)
        mass_p_list.append(mass_p7)
        mass_p_list.append(mass_p8)
        mass_p_list.append(mass_p9)
        mass_p_list.append(mass_p10)
        mass_p_list.append(mass_p11)
        mass_p_list.append(mass_p12)
        mass_p_list.append(mass_p13)
        mass_p_list.append(mass_p14)
        mass_p_list.append(mass_p15)
        print(mass_p_list)

        mass_index_p1 = combobox_p1.current()
        mass_index_p2 = combobox_p2.current()
        mass_index_p3 = combobox_p3.current()
        mass_index_p4 = combobox_p4.current()
        mass_index_p5 = combobox_p5.current()
        mass_index_p6 = combobox_p6.current()
        mass_index_p7 = combobox_p7.current()
        mass_index_p8 = combobox_p8.current()
        mass_index_p9 = combobox_p9.current()
        mass_index_p10 = combobox_p10.current()
        mass_index_p11 = combobox_p011.current()
        mass_index_p12 = combobox_p012.current()
        mass_index_p13 = combobox_p013.current()
        mass_index_p14 = combobox_p014.current()
        mass_index_p15 = combobox_p015.current()
        mass_index_p_list =[]
        mass_index_p_list.append(mass_index_p1)
        mass_index_p_list.append(mass_index_p2)
        mass_index_p_list.append(mass_index_p3)
        mass_index_p_list.append(mass_index_p4)
        mass_index_p_list.append(mass_index_p5)
        mass_index_p_list.append(mass_index_p6)
        mass_index_p_list.append(mass_index_p7)
        mass_index_p_list.append(mass_index_p8)
        mass_index_p_list.append(mass_index_p9)
        mass_index_p_list.append(mass_index_p10)
        mass_index_p_list.append(mass_index_p11)
        mass_index_p_list.append(mass_index_p12)
        mass_index_p_list.append(mass_index_p13)
        mass_index_p_list.append(mass_index_p14)
        mass_index_p_list.append(mass_index_p15)
        print(mass_index_p_list)

        cor_p1 = textBox_p1.get()
        cor_p2 = textBox_p2.get()
        cor_p3 = textBox_p3.get()
        cor_p4 = textBox_p4.get()
        cor_p5 = textBox_p5.get()
        cor_p6 = textBox_p6.get()
        cor_p7 = textBox_p7.get()
        cor_p8 = textBox_p8.get()
        cor_p9 = textBox_p9.get()
        cor_p10 = textBox_p10.get()
        cor_p11 = textBox_p011.get()
        cor_p12 = textBox_p012.get()
        cor_p13 = textBox_p013.get()
        cor_p14 = textBox_p014.get()
        cor_p15 = textBox_p015.get()
        cor_p_list=[]
        cor_p_list.append(cor_p1)
        cor_p_list.append(cor_p2)
        cor_p_list.append(cor_p3)
        cor_p_list.append(cor_p4)
        cor_p_list.append(cor_p5)
        cor_p_list.append(cor_p6)
        cor_p_list.append(cor_p7)
        cor_p_list.append(cor_p8)
        cor_p_list.append(cor_p9)
        cor_p_list.append(cor_p10)
        cor_p_list.append(cor_p11)
        cor_p_list.append(cor_p12)
        cor_p_list.append(cor_p13)
        cor_p_list.append(cor_p14)
        cor_p_list.append(cor_p15)
        print(cor_p_list)
        
        global bln_p1
        global bln_p2
        global bln_p3
        global bln_p4
        global bln_p5
        global bln_p6
        global bln_p7
        global bln_p8
        global bln_p9
        global bln_p10

        size_calc_list = []
        size_error = 0
        if bln_p1.get():
            if ele_list_min_p[0] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[0])
        if bln_p2.get():
            if ele_list_min_p[1] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[1])
        if bln_p3.get():
            if ele_list_min_p[2] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[2])
        if bln_p4.get():
            if ele_list_min_p[3] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[3])
        if bln_p5.get():
            if ele_list_min_p[4] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[4])
        if bln_p6.get():
            if ele_list_min_p[5] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[5])
        if bln_p7.get():
            if ele_list_min_p[6] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[6])
        if bln_p8.get():
            if ele_list_min_p[7] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[7])
        if bln_p9.get():
            if ele_list_min_p[8] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[8])
        if bln_p10.get():
            if ele_list_min_p[9] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[9])
        if bln_p11.get():
            if ele_list_min_p[10] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[10])
        if bln_p12.get():
            if ele_list_min_p[11] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[11])
        if bln_p13.get():
            if ele_list_min_p[12] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[12])
        if bln_p14.get():
            if ele_list_min_p[13] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[13])
        if bln_p15.get():
            if ele_list_min_p[14] == "-":
                size_error = size_error+1
            else:
                size_calc_list.append(ele_list_min_p[14])

        global size_kijun
        size_kijun = combobox_p23.get()

        iso_list_abd_pi = [k for k, x in enumerate(element_list["Atomic Symbol"]) if x == combobox_p23.get()]
        mass_list_abd_pi = []
        for i in range(0,len(iso_list_abd_pi)):
            iso_index_abd_p_i = iso_list_abd_pi[i]
            massnumber = element_list["Mass Number"]
            mass_abd_p_i = massnumber[int(iso_index_abd_p_i)]
            print(mass_abd_p_i)
            mass_list_abd_pi.append(mass_abd_p_i)
        print("masslist")
        print(mass_list_abd_pi)

        for i in range (0,15):
            if combobox_p23.get() == ele_list_min_p[i]:
                cor_pi = cor_p_list[i]
                mass_pi = mass_p_list[i]
        
        global combo23_mass
        combo23_mass = mass_pi

        if ele_p1 != "-" and cor_p1 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p2 != "-" and cor_p2 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p3 != "-" and cor_p3 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p4 != "-" and cor_p4 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p5 != "-" and cor_p5 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p6 != "-" and cor_p6 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p7 != "-" and cor_p7 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p8 != "-" and cor_p8 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p9 != "-" and cor_p9 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p10 != "-" and cor_p10 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p11 != "-" and cor_p11 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p12 != "-" and cor_p12 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p13 != "-" and cor_p13 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p14 != "-" and cor_p14 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif ele_p15 != "-" and cor_p15 == "":
            messagebox.showerror('エラー',"補正項を入力してください")
        elif not size_error == 0:
            messagebox.showerror("エラー","元素のないところにチェックを付けないでください")
        else:
            newwindow_p.destroy()

            newwindow_p2 = tk.Toplevel(root)
            newwindow_p2.geometry("1850x980")
            newwindow_p2.title(u"頂点設定")
            newwindow_p2.configure(bg=color8)

            global data_concat_syuturyoku_hanni_large_p
            data_concat_syuturyoku_hanni_large_p = pd.DataFrame()

            ele_list_p.insert(0,"-")

            frame_p3 = tk.Canvas(newwindow_p2,width=315,height=215,bg=color2)
            frame_p3.place(x=8,y=28)

            frame_p6=tk.Canvas(newwindow_p2,width=71,height=20,bg=color3)
            frame_p6.place(x=16,y=326)
            frame_p7=tk.Canvas(newwindow_p2,width=84,height=20,bg=color4)
            frame_p7.place(x=16,y=406)
            frame_p8=tk.Canvas(newwindow_p2,width=61,height=20,bg=color5)
            frame_p8.place(x=16,y=516)
            frame_p9=tk.Canvas(newwindow_p2,width=84,height=20,bg=color6)
            frame_p9.place(x=16,y=656)

            label_p11 = tk.Label(newwindow_p2,text="A",anchor=tk.CENTER,bg=color2,font=(label_font_jp,10))
            label_p11.place(width=50,height=25,x=12,y=45)
            label_p12 = tk.Label(newwindow_p2,text="B",anchor=tk.CENTER,bg=color2,font=(label_font_jp,10))
            label_p12.place(width=50,height=25,x=12,y=100)
            label_p13 = tk.Label(newwindow_p2,text="C",anchor=tk.CENTER,bg=color2,font=(label_font_jp,10))
            label_p13.place(width=50,height=25,x=12,y=155)

            combobox_p11 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p11.place(width=50,height=25,x=100,y=45)
            combobox_p11.current(0)
            label_p14 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p14.place(width=20,height=25,x=152,y=45)
            combobox_p12 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p12.place(width=50,height=25,x=174,y=45)
            combobox_p12.current(0)
            label_p15 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p15.place(width=20,height=25,x=226,y=45)
            combobox_p13 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p13.place(width=50,height=25,x=248,y=45)
            combobox_p13.current(0)

            combobox_p14 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p14.place(width=50,height=25,x=100,y=100)
            combobox_p14.current(0)
            label_p16 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p16.place(width=20,height=25,x=152,y=100)
            combobox_p15 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p15.place(width=50,height=25,x=174,y=100)
            combobox_p15.current(0)
            label_p17 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p17.place(width=20,height=25,x=226,y=100)
            combobox_p16 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p16.place(width=50,height=25,x=248,y=100)
            combobox_p16.current(0)

            combobox_p17 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p17.place(width=50,height=25,x=100,y=155)
            combobox_p17.current(0)
            label_p18 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p18.place(width=20,height=25,x=152,y=155)
            combobox_p18 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p18.place(width=50,height=25,x=174,y=155)
            combobox_p18.current(0)
            label_p19 = tk.Label(newwindow_p2,text="+",anchor=tk.CENTER,bg=color2,font=(label_font_jp,12))
            label_p19.place(width=20,height=25,x=226,y=155)
            combobox_p19 = ttk.Combobox(newwindow_p2,state='readonly',values=ele_list_min_bar_p)
            combobox_p19.place(width=50,height=25,x=248,y=155)
            combobox_p19.current(0)

            frame_p1=tk.Canvas(newwindow_p2,width=652,height=652,bg=color3)
            frame_p1.place(x=334,y=22)
            frame_p2=tk.Canvas(newwindow_p2,width=610,height=321,bg=color5)
            frame_p2.place(x=996,y=22)
            frame_p4=tk.Canvas(newwindow_p2,width=610,height=321,bg=color6)
            frame_p4.place(x=996,y=353)

            fig_p3 = plt.figure(figsize=(6.42,6.42))
            ax_p3 = fig_p3.add_subplot(111)
            ax_p3.set_aspect('equal', 'datalim')
            plt.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)
            plt.tick_params(bottom=False, left=False, right=False, top=False)
            plt.gca().spines['bottom'].set_visible(False)
            plt.gca().spines['left'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            canvas_p3=FigureCanvasTkAgg(fig_p3,master=newwindow_p2)
            canvas_p3.get_tk_widget().place(x=341,y=29)
            canvas_p3._tkcanvas.place(x=341,y=29)

            fig_p5=plt.figure(figsize=(6,3.11))
            fig_p5.subplots_adjust(bottom=0.2)
            ax_p5=fig_p5.add_subplot(111)
            ax_p5.spines['top'].set_color("w")
            ax_p5.spines['bottom'].set_color("w")
            ax_p5.spines['left'].set_color("w")
            ax_p5.spines['right'].set_color("w")
            ax_p5.tick_params(colors="w")
            canvas_p5=FigureCanvasTkAgg(fig_p5,master=newwindow_p2)
            canvas_p5.get_tk_widget().place(x=1003,y=29)
            canvas_p5._tkcanvas.place(x=1003,y=29)

            fig_p6=plt.figure(figsize=(6,3.11))
            fig_p6.subplots_adjust(bottom=0.2)
            ax_p6=fig_p6.add_subplot(111)
            ax_p6.spines['top'].set_color("w")
            ax_p6.spines['bottom'].set_color("w")
            ax_p6.spines['left'].set_color("w")
            ax_p6.spines['right'].set_color("w")
            ax_p6.tick_params(colors="w")
            canvas_p6=FigureCanvasTkAgg(fig_p6,master=newwindow_p2)
            canvas_p6.get_tk_widget().place(x=1003,y=360)
            canvas_p6._tkcanvas.place(x=1003,y=360)

            fig_p1_5=plt.figure(figsize=(5.4,2.8))
            fig_p1_5.subplots_adjust(bottom=0.2,left=0.15,right=0.75,top=0.9)
            ax_p1_5=fig_p1_5.add_subplot(111)
            ax_p1_5.spines['top'].set_color("w")
            ax_p1_5.spines['bottom'].set_color("w")
            ax_p1_5.spines['left'].set_color("w")
            ax_p1_5.spines['right'].set_color("w")
            ax_p1_5.tick_params(colors="w")
            canvas_p1_5=FigureCanvasTkAgg(fig_p1_5,master=newwindow_p2)
            canvas_p1_5.get_tk_widget().place(x=341,y=691)
            canvas_p1_5._tkcanvas.place(x=341,y=691)

            combobox_p21=ttk.Combobox(newwindow_p2,state='readonly',values=["A","B","C"])
            combobox_p21.place(width=40,height=20,x=110,y=550)
            combobox_p21.current(0)
            combobox_p22=ttk.Combobox(newwindow_p2,state='readonly',values=["A","B","C"])
            combobox_p22.place(width=40,height=20,x=110,y=580)
            combobox_p22.current(0)

            label_p020=tk.Label(newwindow_p2,text="三角グラフ",anchor=tk.W,bg=color3,font=(label_font_jp,10))#,'underline'))
            label_p020.place(width=67,height=16,x=20,y=330)

            label_p021=tk.Label(newwindow_p2,text="三角柱グラフ",anchor=tk.W,bg=color4,font=(label_font_jp,10))#,'underline'))
            label_p021.place(width=80,height=16,x=20,y=410)

            label_p022=tk.Label(newwindow_p2,text="組成分布",anchor=tk.W,bg=color5,font=(label_font_jp,10))#,'underline'))
            label_p022.place(width=57,height=16,x=20,y=520)

            label_p023=tk.Label(newwindow_p2,text="カウント分布",anchor=tk.W,bg=color6,font=(label_font_jp,10))#,'underline'))
            label_p023.place(width=80,height=16,x=20,y=660)

            def click_p2():
                sample_label_long_p = textBox2.get()
                global sample_label_p
                sample_label_p = sample_label_long_p[:-22]

                all_file_p = glob.glob(sample_label_p+"_*_NP_events_large.csv")
                data_num_p = len(all_file_p)

                sample_p = ["{}_{}_NP_events_large".format(sample_label_p,str(i)) for i in range(1,data_num_p+1)]
                print("sample")
                print(sample_p)
                print(len(sample_p))

                global data_p
                global data_use_p
                global data_use_cor_p
                global data_use_size_p
                global data_chouten_zenbu_p
                global data_choutensize_zenbu_p
                global data_choutensize_zentai
                global data_choutensize_kijun

                data_p = pd.DataFrame()
                for run in range (0,len(sample_p)):
                    datasheet_p = "{}.csv".format(sample_p[run])
                    print("\n{}".format(datasheet_p))
                    try:
                        df_p = pd.read_csv(datasheet_p,low_memory=True)
                        data_p = data_p.append(df_p,ignore_index=True)
                    except FileNotFoundError:
                        messagebox.showerror('エラー','選択されたファイルが正しくありません')
                print("data")
                print(data_p)
                global data_concat_syuturyoku_large_p
                data_concat_syuturyoku_large_p = data_p

                data_zero_p = pd.DataFrame(np.zeros_like(data_p["t_elapsed_Buf"]))

                data_use_p = data_p["'[{}{}]+'".format(mass_p1,ele_p1)]                
                for i in range(1,len(ele_list_min_ini_p)):
                    data_i_p = data_p["'[{}{}]+'".format(mass_p_list[i],ele_p_list[i])]
                    data_use_p = pd.concat([data_use_p,data_i_p],axis=1)
                data_use_p = pd.concat([data_use_p,data_zero_p],axis=1)


                data_use_cor_p = data_p["'[{}{}]+'".format(mass_p1,ele_p1)]/float(cor_p1)
                for i in range(1,len(ele_list_min_ini_p)):
                    data_i_cor_p = data_p["'[{}{}]+'".format(mass_p_list[i],ele_p_list[i])]/float(cor_p_list[i])
                    data_use_cor_p = pd.concat([data_use_cor_p,data_i_cor_p],axis=1)
                data_use_cor_p = pd.concat([data_use_cor_p,data_zero_p],axis=1)

                Iso_Compo_list = element_list["Isotopic Composition"]
                list_number_kijun = mass_list_abd_pi.index(int(combo23_mass))
                index_pi = iso_list_abd_pi[int(list_number_kijun)]
                iso_compo_pi = Iso_Compo_list[int(index_pi)]
                print("iso_compo_pi")
                print(iso_compo_pi)
                
                data_use_size_p = data_p["'[{}{}]+'".format(mass_p1,ele_p1)]/float(iso_compo_pi)*float(cor_pi)/float(cor_p1)           
                for i in range(1,len(ele_list_min_ini_p)):
                    data_i_size_p = data_p["'[{}{}]+'".format(mass_p_list[i],ele_p_list[i])]/float(iso_compo_pi)*float(cor_pi)/float(cor_p_list[i])
                    data_use_size_p = pd.concat([data_use_size_p,data_i_size_p],axis=1)
                data_use_size_p = pd.concat([data_use_size_p,data_zero_p],axis=1)

                data_use_p.columns = ele_list_min_zero_p
                data_use_cor_p.columns = ele_list_min_zero_p
                data_use_size_p.columns = ele_list_min_zero_p

                data_chouten1_p = data_use_cor_p[combobox_p11.get()]+data_use_cor_p[combobox_p12.get()]+data_use_cor_p[combobox_p13.get()]
                data_chouten2_p = data_use_cor_p[combobox_p14.get()]+data_use_cor_p[combobox_p15.get()]+data_use_cor_p[combobox_p16.get()]
                data_chouten3_p = data_use_cor_p[combobox_p17.get()]+data_use_cor_p[combobox_p18.get()]+data_use_cor_p[combobox_p19.get()]
                data_chouten_zenbu_p = pd.concat([data_chouten1_p,data_chouten2_p,data_chouten3_p],axis=1)
                data_chouten_zenbu_p.columns = ["CT1","CT2","CT3"]

                data_choutensize1_p = data_use_size_p[combobox_p11.get()]+data_use_size_p[combobox_p12.get()]+data_use_size_p[combobox_p13.get()]
                data_choutensize2_p = data_use_size_p[combobox_p14.get()]+data_use_size_p[combobox_p15.get()]+data_use_size_p[combobox_p16.get()]
                data_choutensize3_p = data_use_size_p[combobox_p17.get()]+data_use_size_p[combobox_p18.get()]+data_use_size_p[combobox_p19.get()]
                data_choutensize_zenbu_p = pd.concat([data_choutensize1_p,data_choutensize2_p,data_choutensize3_p],axis=1)
                data_choutensize_zenbu_p.columns = ["CT1","CT2","CT3"]

                data_choutensize_zentai = 0
                for i in range(0,len(size_calc_list)):
                    data_choutensize_zentai = data_choutensize_zentai + data_use_size_p[size_calc_list[i]]

                data_choutensize_kijun = data_use_size_p[size_kijun]

            def click_pC():

                plt.rcParams["font.family"] = graph_font

                SUM = data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                per = SUM/data_choutensize_zentai*100
                data_concat_p = pd.concat([pX,pY,pZ,SUM,per],axis=1)
                data_concat_p.columns = ['X','Y','Z','SUM','per']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','X']

                zmax_p = data_z_p['SUM'].max()

                amari_p = zmax_p%10

                if amari_p == 0:
                    baisuu_p = zmax_p
                else:
                    baisuu_p = zmax_p-amari_p+10

                if textBox_p11.get() == "":
                    tB_p11 = "0"
                else:
                    tB_p11 = textBox_p11.get()

                if textBox_p12.get() == "":
                    tB_p12 = baisuu_p
                else:
                    tB_p12 = textBox_p12.get()

                if tB_p11 == "0" and tB_p12 == baisuu_p:
                    data_concat_0_tri_p = data_concat_p
                    data_gattai_0_p = data_chouten_zenbu_p
                    data_concat_syuturyoku_z_large_p = data_concat_syuturyoku_large_p
                else:
                    data_concat_0_tri_p = data_concat_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]
                    data_gattai_0_p = data_chouten_zenbu_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]
                    data_concat_syuturyoku_z_large_p = data_concat_syuturyoku_large_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]
                
                if textBox_p17.get()=="":
                    tB_p17="0"
                else:
                    tB_p17=textBox_p17.get()

                if textBox_p18.get()=="":
                    tB_p18="100"
                else:
                    tB_p18=textBox_p18.get()

                if textBox_p19.get()=="":
                    tB_p19="0"
                else:
                    tB_p19=textBox_p19.get()

                if textBox_p20.get()=="":
                    tB_p20="100"
                else:
                    tB_p20=textBox_p20.get()
                    
                if textBox_p21.get()=="":
                    tB_p21="0"
                else:
                    tB_p21=textBox_p21.get()

                if textBox_p22.get()=="":
                    tB_p22="100"
                else:
                    tB_p22=textBox_p22.get()

                global data_concat_syuturyoku_hanni_large_p
                data_concat0_tri_hanni_p = data_concat_0_tri_p[(data_concat_0_tri_p["X"]>=int(tB_p17))&(data_concat_0_tri_p["X"]<=int(tB_p18))&(data_concat_0_tri_p["Y"]>=int(tB_p19))&(data_concat_0_tri_p["Y"]<=int(tB_p20))&(data_concat_0_tri_p["Z"]>=int(tB_p21))&(data_concat_0_tri_p["Z"]<=int(tB_p22))]
                data_gattai_hanni_p = data_gattai_0_p[(data_concat_0_tri_p["X"]>=int(tB_p17))&(data_concat_0_tri_p["X"]<=int(tB_p18))&(data_concat_0_tri_p["Y"]>=int(tB_p19))&(data_concat_0_tri_p["Y"]<=int(tB_p20))&(data_concat_0_tri_p["Z"]>=int(tB_p21))&(data_concat_0_tri_p["Z"]<=int(tB_p22))]
                data_concat_syuturyoku_hanni_large_p = data_concat_syuturyoku_z_large_p[(data_concat_0_tri_p["X"]>=int(tB_p17))&(data_concat_0_tri_p["X"]<=int(tB_p18))&(data_concat_0_tri_p["Y"]>=int(tB_p19))&(data_concat_0_tri_p["Y"]<=int(tB_p20))&(data_concat_0_tri_p["Z"]>=int(tB_p21))&(data_concat_0_tri_p["Z"]<=int(tB_p22))]

                mean_data_p = data_concat0_tri_hanni_p.mean()
                print("平均")
                print(mean_data_p)

                std_data_p = np.std(data_concat0_tri_hanni_p, ddof=1) #標本標準偏差
                print("標準偏差")
                print(std_data_p)

                ste_data_p = std_data_p/(len(data_concat0_tri_hanni_p)**(1/2))
                print("標準誤差")
                print(ste_data_p)
                
                h = np.sqrt(3.0)*0.5

                global data_plot_p3
                plotx_p3 = (100-data_concat0_tri_hanni_p['Y'])/100-data_concat0_tri_hanni_p['X']/200
                ploty_p3 = h*data_concat0_tri_hanni_p['X']/100
                data_plot_p3 = pd.concat([plotx_p3,ploty_p3],axis=1)
                data_plot_p3.columns = ['x3','y3']

                global data_click_p
                data_click_p = pd.concat([data_concat0_tri_hanni_p,data_plot_p3],axis=1)
                data_click_p.columns=["a'","b'","c'","sum","per","x","y"]

                ax_p3.cla()
                fig_p3.set_facecolor(color_fig)
                ax_p3.set_facecolor(color_fig)

                for i in range(1,10):
                    ax_p3.plot([i/20.0, 1.0-i/20.0],[h*i/10.0, h*i/10.0], linestyle='dashed',color='gray', lw=0.5,zorder=2)
                    ax_p3.plot([i/20.0, i/10.0],[h*i/10.0, 0.0], linestyle='dashed',color='gray', lw=0.5,zorder=2)
                    ax_p3.plot([0.5+i/20.0, i/10.0],[h*(1.0-i/10.0), 0.0], linestyle='dashed',color='gray', lw=0.5,zorder=2)
                
                ax_p3.plot([0.0, 1.0],[0.0, 0.0], color=color_axes, lw=2,zorder=2)
                ax_p3.plot([0.0, 0.5],[0.0, h], color=color_axes, lw=2,zorder=2)
                ax_p3.plot([1.0, 0.5],[0.0, h], color=color_axes, lw=2,zorder=2)   

                chouten1_text = combobox_p11.get()
                chouten2_text = combobox_p14.get()
                chouten3_text = combobox_p17.get()

                ax_p3.text(0.455, h+0.0283, chouten1_text, fontsize=22, color=color_axes, va="center")
                ax_p3.text(-0.1, -0.02, chouten2_text, fontsize=22, color=color_axes, va="center")
                ax_p3.text(1.02, -0.02, chouten3_text, fontsize=22, color=color_axes, va="center")

                for i in range(1,10):
                    ax_p3.text(0.5+(10-i)/20.0+0.016, h*(1.0-(10-i)/10.0), '%d0' % i, fontsize=17, color=color_axes)
                    ax_p3.text((10-i)/20.0-0.082, h*(10-i)/10.0, '%d0' % i, fontsize=17, color=color_axes)
                    ax_p3.text(i/10.0-0.03, -0.06, '%d0' % i, fontsize=17, color=color_axes)

                ax_p3.text(-0.15,1,"Number of particles:"+str(len(data_plot_p3.dropna())),fontsize=14, color=color_axes)

                if tri_plot_color == "ON":
                    hi_list_p = []
                    for i in range(0,10):
                        for j in range(0,10-i):
                            data_concat_hi_p1 = data_concat0_tri_hanni_p[(data_concat0_tri_hanni_p["X"]>=i*10)&(data_concat0_tri_hanni_p["X"]<=(i+1)*10)&(data_concat0_tri_hanni_p["Y"]>=(9-i-j)*10)&(data_concat0_tri_hanni_p["Y"]<=(10-i-j)*10)&(data_concat0_tri_hanni_p["Z"]>=j*10)&(data_concat0_tri_hanni_p["Z"]<=(j+1)*10)]
                            hi_p1=len(data_concat_hi_p1)
                            hi_list_p.append(hi_p1)
                    for i in range(0,10):
                        for j in range(0,10-i-1):
                            data_concat_hi_p2 = data_concat0_tri_hanni_p[(data_concat0_tri_hanni_p["X"]>=i*10)&(data_concat0_tri_hanni_p["X"]<=(i+1)*10)&(data_concat0_tri_hanni_p["Y"]>=(8-i-j)*10)&(data_concat0_tri_hanni_p["Y"]<=(9-i-j)*10)&(data_concat0_tri_hanni_p["Z"]>=j*10)&(data_concat0_tri_hanni_p["Z"]<=(j+1)*10)]
                            hi_p2=len(data_concat_hi_p2)
                            hi_list_p.append(hi_p2)

                    for i in range(0,10):
                        for j in range(0,10-i):
                            data_concat_hi_p = data_concat0_tri_hanni_p[(data_concat0_tri_hanni_p["X"]>=i*10)&(data_concat0_tri_hanni_p["X"]<=(i+1)*10)&(data_concat0_tri_hanni_p["Y"]>=(9-i-j)*10)&(data_concat0_tri_hanni_p["Y"]<=(10-i-j)*10)&(data_concat0_tri_hanni_p["Z"]>=j*10)&(data_concat0_tri_hanni_p["Z"]<=(j+1)*10)]
                            hi_p=len(data_concat_hi_p)
                            if hi_p >= hi_bord_1:
                                hi_color_p = "red"
                            elif hi_p >= hi_bord_2:
                                hi_color_p = "darkorange"
                            elif hi_p >= hi_bord_3:
                                hi_color_p = "lawngreen"
                            elif hi_p >= hi_bord_4:
                                hi_color_p = "green"
                            elif hi_p >= hi_bord_5:
                                hi_color_p = "blue"
                            elif hi_p > 0:
                                hi_color_p = "lightskyblue"
                            else:
                                hi_color_p = "#ffffff00"
                            x_1 = i/20+j/10
                            x_2 = i/20+(j+1)/10
                            x_3 = i/20+j/10+1/20
                            y_1 = i*h/10
                            y_2 = i*h/10
                            y_3 = i*h/10+h/10
                            poly = pat.Polygon(xy = [(x_1,y_1),(x_2,y_2),(x_3,y_3)],fc=hi_color_p,ec=None,alpha=0.5,zorder=1)
                            ax_p3.add_patch(poly)

                    for i in range(0,10):
                        for j in range(0,10-i-1):
                            data_concat_hi_p = data_concat0_tri_hanni_p[(data_concat0_tri_hanni_p["X"]>=i*10)&(data_concat0_tri_hanni_p["X"]<=(i+1)*10)&(data_concat0_tri_hanni_p["Y"]>=(8-i-j)*10)&(data_concat0_tri_hanni_p["Y"]<=(9-i-j)*10)&(data_concat0_tri_hanni_p["Z"]>=j*10)&(data_concat0_tri_hanni_p["Z"]<=(j+1)*10)]
                            hi_p=len(data_concat_hi_p)
                            if hi_p >= hi_bord_1: #hi/hi_max >= 0.8:
                                hi_color_p = "red"
                            elif hi_p >= hi_bord_2: #hi/hi_max >= 0.6:
                                hi_color_p = "darkorange"
                            elif hi_p >= hi_bord_3: #hi/hi_max >= 0.4:
                                hi_color_p = "lawngreen"
                            elif hi_p >= hi_bord_4: #hi/hi_max >= 0.2:
                                hi_color_p = "green"
                            elif hi_p >= hi_bord_5: #hi/hi_max >= 0.1:
                                hi_color_p = "blue"
                            elif hi_p > 0: #hi/hi_max > 0:
                                hi_color_p = "lightskyblue"
                            else:
                                hi_color_p = "#ffffff00"
                            x_1 = i/20+j/10+3/20
                            x_2 = i/20+(j+1)/10
                            x_3 = i/20+j/10+1/20
                            y_1 = i*h/10+h/10
                            y_2 = i*h/10
                            y_3 = i*h/10+h/10
                            poly = pat.Polygon(xy = [(x_1,y_1),(x_2,y_2),(x_3,y_3)],fc=hi_color_p,ec=None,alpha=0.5,zorder=1)
                            ax_p3.add_patch(poly)
                elif tri_plot_color == "OFF":
                    print("OFF")
                else:
                    print("ERROR")
                
                ax_p3.scatter(data_plot_p3["x3"],data_plot_p3["y3"],c=color_plot,alpha=alpha_plot,s=size_plot,zorder=3)

                canvas_p3.draw()                     

            def click_pA():

                plt.rcParams["font.family"] = graph_font

                SUM = data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                data_concat_p = pd.concat([pX,pY,pZ,SUM],axis=1)
                data_concat_p.columns = ['X','Y','Z','SUM']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','X']

                h = np.sqrt(3.0)*0.5

                zmax_p = data_z_p['SUM'].max()

                amari_p = zmax_p%10
                
                if amari_p == 0:
                    baisuu_p = zmax_p
                else:
                    baisuu_p = zmax_p-amari_p+10

                NP1_p = np.linspace(0,baisuu_p,6,dtype=int)

                plotx_p1 = (100-data_concat_p['Y'])/100-data_concat_p['X']/200
                ploty_p1 = h*data_concat_p['X']/100
                plotz_p1 = data_z_p['SUM']/baisuu_p
                data_plot_p1 = pd.concat([plotx_p1,ploty_p1,plotz_p1],axis=1)
                data_plot_p1.columns = ['x1','y1','z1']

                fig_p1 = plt.figure(figsize=(7.5,7.5))
                ax_p1 = fig_p1.add_subplot(111,projection='3d')
                ax_p1.set_xticks([])
                ax_p1.set_yticks([])
                ax_p1.set_zticks([])
                plt.axis('off')
                ax_p1.view_init(elev=12,azim=-69)

                fig_p1.set_facecolor(color_fig)
                ax_p1.set_facecolor(color_fig)

                for i in range(1,5):
                    ax_p1.plot([i*2/20.0, 1.0-i*2/20.0],[h*2*i/10.0, h*i*2/10.0],[0,0],color='gray', lw=0.5)
                    ax_p1.plot([i*2/20.0, i*2/10.0],[h*i*2/10.0, 0.0],[0,0], color='gray', lw=0.5)
                    ax_p1.plot([0.5+i*2/20.0, i*2/10.0],[h*(1.0-i*2/10.0), 0.0],[0,0],color='gray', lw=0.5)

                ax_p1.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
                ax_p1.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p1.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p1.plot([0.0, 0.0],[0.0, 0.0],[0,1], color=color_axes, lw=2)
                ax_p1.plot([1, 1],[0,0],[0,1], color=color_axes, lw=2)
                ax_p1.plot([0.5, 0.5],[h, h],[0,1], color=color_axes, lw=2)
                ax_p1.text(0.5, h+0.15,-0.1, combobox_p11.get(), fontsize=20,ha="center", color=color_axes)
                ax_p1.text(-0.15*h, -0.15/2,-0.1, combobox_p14.get(), fontsize=20,ha="center", color=color_axes)
                ax_p1.text(1+0.15*h, -0.15/2,-0.1, combobox_p17.get(), fontsize=20,ha="center", color=color_axes)

                for i in range(1,5):
                    ax_p1.plot([2*i/20.0, 1.0-2*i/20.0],[2*h*i/10.0, 2*h*i/10.0],[1,1],color='gray', lw=0.5)
                    ax_p1.plot([2*i/20.0, 2*i/10.0],[2*h*i/10.0, 0.0],[1,1], color='gray', lw=0.5)
                    ax_p1.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[1,1], color='gray', lw=0.5)
    
                ax_p1.plot([0.0, 1.0],[0.0, 0.0],[1,1], color=color_axes, lw=2)
                ax_p1.plot([0.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)
                ax_p1.plot([1.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)

                ax_p1.plot([0,0],[0,0],[0,1], color=color_axes, lw=2)

                for i in range(0,6):
                   ax_p1.plot([0,-0.02],[0,-0.02],[i/5,i/5], color=color_axes, lw=2)

                for i in range(0,6):
                   ax_p1.text(-0.078,-0.078,i/5-0.01,NP1_p[i],fontsize=20,ha="right", color=color_axes)

                if size_hosei_p == "":
                    ax_p1.text(-0.08,-0.08,1.2,"Counts",ha="right",fontsize=20,color=color_axes) 
                else:
                    ax_p1.text(-0.08,-0.08,1.2,"Size (nm)",ha="right",fontsize=20,color=color_axes)
    
                ax_p1.set_xlim(0,1)
                ax_p1.set_ylim(0,1)
                ax_p1.set_zlim(0,1)

                #プロット
                ax_p1.scatter(data_plot_p1["x1"],data_plot_p1["y1"],data_plot_p1["z1"],c=color_plot,alpha=alpha_plot,depthshade=False,s=size_plot)

                #範囲指定のところに値が入っていた場合に全体表示のどの範囲にあたるかを示す
                if not textBox_p13.get() == "" or not textBox_p14.get() == "" :
                    if textBox_p13.get()=="":
                       tB_p13="0"
                    else:
                       tB_p13=textBox_p13.get()

                    if textBox_p14.get()=="":
                       tB_p14=baisuu_p
                    else:
                       tB_p14=textBox_p14.get()

                    tB_p13_range = int(tB_p13)/baisuu_p
                    tB_p14_range = int(tB_p14)/baisuu_p   

                    ax_p1.plot([0.0, 1.0],[0.0, 0.0],[tB_p13_range,tB_p13_range], color="red", lw=1.5)
                    ax_p1.plot([0.0, 0.5],[0.0, h],[tB_p13_range,tB_p13_range], color="red", lw=1.5)
                    ax_p1.plot([1.0, 0.5],[0.0, h],[tB_p13_range,tB_p13_range], color="red", lw=1.5)
                    ax_p1.plot([0.0, 0.0],[0.0, 0.0],[tB_p13_range,tB_p14_range], color="red", lw=1.5)
                    ax_p1.plot([1, 1],[0,0],[tB_p13_range,tB_p14_range], color="red", lw=1.5)
                    ax_p1.plot([0.5, 0.5],[h, h],[tB_p13_range,tB_p14_range], color="red", lw=1.5)
                    ax_p1.plot([0.0, 1.0],[0.0, 0.0],[tB_p14_range,tB_p14_range], color="red", lw=1.5)
                    ax_p1.plot([0.0, 0.5],[0.0, h],[tB_p14_range,tB_p14_range], color="red", lw=1.5)
                    ax_p1.plot([1.0, 0.5],[0.0, h],[tB_p14_range,tB_p14_range], color="red", lw=1.5)        

                fig_p1.show()

            def click_pD(): 

                plt.rcParams["font.family"] = graph_font

                SUM = data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                data_concat_p = pd.concat([pX,pY,pZ,SUM],axis=1)
                data_concat_p.columns = ['X','Y','Z','SUM']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','X']

                h = np.sqrt(3.0)*0.5
     
                #座標決め
                zmax_p = data_z_p['SUM'].max()

                amari_p = zmax_p%10
                
                if amari_p == 0:
                    baisuu_p = zmax_p
                else:
                    baisuu_p = zmax_p-amari_p+10

                plotx_p4=(100-data_concat_p['Y'])/100-data_concat_p['X']/200
                ploty_p4=h*data_concat_p['X']/100

                if textBox_p13.get()=="":
                   tB_p13="0"
                else:
                   tB_p13=textBox_p13.get()

                if textBox_p14.get()=="":
                   tB_p14=baisuu_p
                else:
                   tB_p14=textBox_p14.get()

                if tB_p13=="0" and tB_p14==baisuu_p:
                   zvalue_p=baisuu_p
                   plotz_p4=data_z_p['SUM']/zvalue_p
                   NP4_p=np.linspace(0,baisuu_p,6,dtype=int)
                else:
                   zvalue_p=int(tB_p14)-int(tB_p13)
                   plotz_p4=(data_z_p['SUM']-int(tB_p13))/zvalue_p
                   NP4_p=np.linspace(int(tB_p13),int(tB_p14),6,dtype=int)

                data_plot_p4=pd.concat([plotx_p4,ploty_p4,plotz_p4],axis=1)
                data_plot_p4.columns=['x4','y4','z4']
                data_plot_p40=data_plot_p4[(data_plot_p4["z4"]>=0)&(data_plot_p4["z4"]<=1)]

                #プロット場所決め
                fig_p4=plt.figure(figsize=(7.5,7.5))
                ax_p4=fig_p4.add_subplot(111,projection='3d')
                ax_p4.set_xticks([])
                ax_p4.set_yticks([])
                ax_p4.set_zticks([])
                plt.axis('off')
                ax_p4.view_init(elev=12,azim=-69)

                #三角グラフ描き
                fig_p4.set_facecolor(color_fig)
                ax_p4.set_facecolor(color_fig)

                for i in range(1,5):
                    ax_p4.plot([2*i/20.0, 1.0-2*i/20.0],[h*2*i/10.0, h*2*i/10.0],[0,0],color='gray', lw=0.5)
                    ax_p4.plot([2*i/20.0, 2*i/10.0],[h*2*i/10.0, 0.0],[0,0], color='gray', lw=0.5)
                    ax_p4.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[0,0], color='gray', lw=0.5)

                ax_p4.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
                ax_p4.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p4.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p4.plot([0.0, 0.0],[0.0, 0.0],[0,1], color=color_axes, lw=2)
                ax_p4.plot([1, 1],[0,0],[0,1], color=color_axes, lw=2)
                ax_p4.plot([0.5, 0.5],[h, h],[0,1], color=color_axes, lw=2)
                ax_p4.text(0.5, h+0.15,-0.1, combobox_p11.get(), fontsize=20,ha="center", color=color_axes)
                ax_p4.text(-0.15*h, -0.15/2,-0.1, combobox_p14.get(), fontsize=20,ha="center", color=color_axes)
                ax_p4.text(1+0.15*h, -0.15/2,-0.1, combobox_p17.get(), fontsize=20,ha="center", color=color_axes)

                for i in range(1,5):
                    ax_p4.plot([2*i/20.0, 1.0-2*i/20.0],[2*h*i/10.0, 2*h*i/10.0],[1,1],color='gray', lw=0.5)
                    ax_p4.plot([2*i/20.0, 2*i/10.0],[2*h*i/10.0, 0.0],[1,1], color='gray', lw=0.5)
                    ax_p4.plot([0.5+2*i/20.0, 2*i/10.0],[h*(1.0-2*i/10.0), 0.0],[1,1], color='gray', lw=0.5)

                ax_p4.plot([0.0, 1.0],[0.0, 0.0],[1,1], color=color_axes, lw=2)
                ax_p4.plot([0.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)
                ax_p4.plot([1.0, 0.5],[0.0, h],[1,1], color=color_axes, lw=2)

                ax_p4.plot([0,0],[0,0],[0,1], color=color_axes, lw=2)

                #for i in range(0,6):
                #   ax_p4.plot([0,-0.02],[0,-0.02],[i/5,i/5], color=color_axes, lw=2)

                #for i in range(0,6):
                #   ax_p4.text(-0.078,-0.078,i/5-0.01,NP4_p[i],fontsize=20,ha="right", color=color_axes)

                if size_hosei_p == "":
                    ax_p4.text(-0.08,-0.08,1.2,"Counts",ha="right",fontsize=20,color=color_axes) 
                else:
                    ax_p4.text(-0.08,-0.08,1.2,"Size (nm)",ha="right",fontsize=20,color=color_axes)

                ax_p4.set_xlim(0,1)
                ax_p4.set_ylim(0,1)
                ax_p4.set_zlim(0,1)

                #プロット
                ax_p4.scatter(data_plot_p40["x4"],data_plot_p40["y4"],data_plot_p40["z4"],c=color_plot,alpha=alpha_plot,depthshade=False,s=size_plot)

                fig_p4.show()

            def click_pE():

                plt.rcParams["font.family"] = graph_font

                SUM = data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                data_concat_p = pd.concat([pX,pY,pZ,SUM],axis=1)
                data_concat_p.columns = ['X','Y','Z','SUM']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','X']

                zmax_p = data_z_p['SUM'].max()

                amari_p = zmax_p%10

                if amari_p == 0:
                    baisuu_p = zmax_p
                else:
                    baisuu_p = zmax_p-amari_p+10

                if textBox_p11.get() == "":
                    tB_p11 = "0"
                else:
                    tB_p11 = textBox_p11.get()

                if textBox_p12.get() == "":
                    tB_p12 = baisuu_p
                else:
                    tB_p12 = textBox_p12.get()

                if tB_p11 == "0" and tB_p12 == baisuu_p:
                    data_concat0_p = data_concat_p
                else:
                    data_concat0_p = data_concat_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]

                if textBox_p15.get()=="":
                    tB_p15="0"
                else:
                    tB_p15=textBox_p15.get()

                if textBox_p16.get()=="":
                    tB_p16="100"
                else:
                    tB_p16=textBox_p16.get()

                if combobox_p21.get()=="A":
                    if tB_p15=="0" and tB_p16=="100":
                        data_concat00_p=data_concat0_p
                    else:
                        data_concat00_p=data_concat0_p[(data_concat0_p['X']>=int(tB_p15))&(data_concat0_p['X']<=int(tB_p16))]

                    if combobox_p22.get()=="B":
                        count_p=data_concat00_p['Y']/(data_concat00_p['Y']+data_concat00_p['Z'])*100
                        perelement_p="B"
                    elif combobox_p22.get()=="C":
                        count_p=data_concat00_p['Z']/(data_concat00_p['Y']+data_concat00_p['Z'])*100
                        perelement_p="C"   
                    elif combobox_p22.get()=="":
                        messagebox.showerror('エラー','割合を見る元素を指定してください') 
                    else:
                        messagebox.showerror('エラー','重複せずに選択してください') 
                elif combobox_p21.get()=="B":        
                    if tB_p15=="0" and tB_p16=="100":
                        data_concat00_p=data_concat0_p
                    else:
                        data_concat00_p=data_concat0_p[(data_concat0_p['Y']>=int(tB_p15))&(data_concat0_p['Y']<=int(tB_p16))]
      
                    if combobox_p22.get()=="C":
                        count_p=data_concat00_p['Z']/(data_concat00_p['Z']+data_concat00_p['X'])*100
                        perelement_p="C"
                    elif combobox_p22.get()=="A":
                        count_p=data_concat00_p['X']/(data_concat00_p['Z']+data_concat00_p['X'])*100
                        perelement_p="A"
                    elif combobox_p22.get()=="":
                        messagebox.showerror('エラー','割合を見る元素を指定してください')       
                    else:
                        messagebox.showerror('エラー','重複せずに選択してください') 
                elif combobox_p21.get()=="C":
                    if tB_p15=="0" and tB_p16=="100":
                        data_concat00_p=data_concat0_p
                    else:
                        data_concat00_p=data_concat0_p[(data_concat0_p['Z']>=int(tB_p15))&(data_concat0_p['Z']<=int(tB_p16))]
      
                    if combobox_p22.get()=="A":
                        count_p=data_concat00_p['X']/(data_concat00_p['X']+data_concat00_p['Y'])*100
                        perelement_p="A"
                    elif combobox_p22.get()=="B":
                        count_p=data_concat00_p['Y']/(data_concat00_p['X']+data_concat00_p['Y'])*100 
                        perelement_p="A"
                    elif combobox_p22.get()=="":
                        messagebox.showerror('error','割合を見る元素を指定してください') 
                    else:
                        messagebox.showerror('error','重複せずに選択してください') 
                elif combobox_p21.get()=="":
                    messagebox.showerror('error','範囲を固定する元素を指定してください')

                data_z_c_p=pd.concat([count_p,data_z_SUM_p],axis=1)
                data_z_c_p.columns=['X','SUM']     
                data_z_0_p = data_z_c_p.dropna(how="any")

                tB_p333=combobox_p20.get()

                times_p=100/int(tB_p333)

                Num4_X_p=[]
                for i in range(0,int(times_p)):
                    Num4_X_p.append(np.count_nonzero((count_p>=i*int(tB_p333))&(count_p<=(i+1)*int(tB_p333))))

                ax_p5.cla()

                plt.rcParams["font.family"] = graph_font

                fig_p5.set_facecolor(color_fig)
                ax_p5.set_facecolor(color_fig)

                ax_p5.spines['top'].set_color(color_axes)
                ax_p5.spines['bottom'].set_color(color_axes)
                ax_p5.spines['left'].set_color(color_axes)
                ax_p5.spines['right'].set_color(color_axes)
                ax_p5.tick_params(colors=color_axes)

                left=np.linspace(0,100-int(tB_p333),int(times_p),dtype=int)
                height=np.array(Num4_X_p)
                if bar_style == "枠あり":   
                    ax_p5.bar(left,height,width=int(tB_p333),color=color_bar,linewidth=1,edgecolor=color_axes,align="edge",zorder=2)
                elif bar_style == "枠なし":
                    ax_p5.bar(left,height,width=int(tB_p333)* 0.8,color=color_bar,align="edge",zorder=2)
                ax_p5.grid(b=True,which='major',axis='y',color=color_axes,linewidth=0.5,zorder=1)
                ax_p5.set_xlabel("%"+perelement_p,fontname=graph_font,color=color_axes,weight=weight_font)
                ax_p5.set_ylabel("Number of particles",fontname=graph_font,color=color_axes,weight=weight_font)

                canvas_p5.draw()

            def click_pF():
                plt.rcParams["font.family"] = graph_font
            
                SUM = data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                data_concat_p = pd.concat([pX,pY,pZ],axis=1)
                data_concat_p.columns = ['X','Y','Z']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai #いらない列
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai #いらない列
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai#いらない列
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','x']

                data_concat_z_p=pd.concat([data_concat_p,data_z_p],axis=1)

                if textBox_p17.get()=="":
                    tB_p17="0"
                else:
                    tB_p17=textBox_p17.get()
            
                if textBox_p18.get()=="":
                    tB_p18="100"
                else:
                    tB_p18=textBox_p18.get()

                if textBox_p19.get()=="":
                    tB_p19="0"
                else:
                    tB_p19=textBox_p19.get()
        
                if textBox_p20.get()=="":
                    tB_p20="100"
                else:
                    tB_p20=textBox_p20.get()

                if textBox_p21.get()=="":
                    tB_p21="0"
                else:
                    tB_p21=textBox_p21.get()

                if textBox_p22.get()=="":
                    tB_p22="100"
                else:
                    tB_p22=textBox_p22.get()

                data_concat_z_p_drop = data_concat_z_p.dropna(how="any")
                data_concat0_p=data_concat_z_p_drop[(data_concat_z_p_drop["X"]>=int(tB_p17))&(data_concat_z_p_drop["X"]<=int(tB_p18))&(data_concat_z_p_drop["Y"]>=int(tB_p19))&(data_concat_z_p_drop["Y"]<=int(tB_p20))&(data_concat_z_p_drop["Z"]>=int(tB_p21))&(data_concat_z_p_drop["Z"]<=int(tB_p22))]
                zmax0_p = data_concat0_p['SUM'].max()
                amari_p=zmax0_p%10
                if amari_p==0:
                   baisuu_p=zmax0_p
                else:
                   baisuu_p=zmax0_p-amari_p+10

                if textBox_p23.get()=="":
                    tB_p23="0"
                else:
                    tB_p23=textBox_p23.get()

                if textBox_p24.get()=="":
                    tB_p24=baisuu_p
                else:
                    tB_p24=textBox_p24.get()
            
                if textBox_p25.get()=="":
                    tB_p25=5
                else:
                    tB_p25=textBox_p25.get()

                if len(data_concat0_p)==0:
                    if int(tB_p21)<100-(int(tB_p18)+int(tB_p20)) or int(tB_p22)>100-(int(tB_p17)+int(tB_p19)):
                        messagebox.showerror('エラー','範囲が間違っています')
                    else:
                        messagebox.showerror('エラー','範囲内にデータがありません')
                else:
                    pltmin_p=int(tB_p23)
                    amari0_p=int(tB_p24)%10
                    if amari0_p==0:
                        pltmax_p=int(tB_p24)
                    else:
                        pltmax_p=int(tB_p24)-amari0_p+10

                    ax_p6.cla()

                    fig_p6.set_facecolor(color_fig)
                    ax_p6.set_facecolor(color_fig)

                    ax_p6.spines['top'].set_color(color_axes)
                    ax_p6.spines['bottom'].set_color(color_axes)
                    ax_p6.spines['left'].set_color(color_axes)
                    ax_p6.spines['right'].set_color(color_axes)
                    ax_p6.tick_params(colors=color_axes)
    
                    if bar_style == "枠あり":
                        ax_p6.hist(data_concat0_p["SUM"],bins=np.arange(pltmin_p,pltmax_p+int(tB_p25),int(tB_p25)),color=color_bar,linewidth=1,edgecolor=color_axes,zorder=2)
                    elif bar_style == "枠なし":
                        ax_p6.hist(data_concat0_p["SUM"],bins=np.arange(pltmin_p,pltmax_p+int(tB_p25),int(tB_p25)),width = int(tB_p25) * 0.8,color=color_bar,zorder=2)
                    ax_p6.grid(b=True,which='major',axis='y',color=color_axes,linewidth=0.5,zorder=1)
                    
                    if size_hosei_p == "":
                        ax_p6.set_xlabel("Counts",fontname=graph_font,color=color_axes,weight=weight_font)
                    else:
                        ax_p6.set_xlabel("Size (nm)",fontname=graph_font,color=color_axes,weight=weight_font)
                    ax_p6.set_ylabel("Number of particles",fontname=graph_font,color=color_axes,weight=weight_font)

                    data_concat0_small = data_concat0_p[(data_concat0_p["SUM"]>=int(tB_p23))&(data_concat0_p["SUM"]<=int(tB_p24))]

                    ax_p6_pos = ax_p6.get_position()
                    fig_p6.text(ax_p6_pos.x1-0.24,ax_p6_pos.y1+0.04,"{} particlespar".format(len(data_concat0_p)),fontsize=14, color=color_fig,backgroundcolor=color_fig)
                    fig_p6.text(ax_p6_pos.x1-0.24,ax_p6_pos.y1+0.04,"{} particles".format(len(data_concat0_small)),fontsize=14, color=color_axes,backgroundcolor=color_fig)

                    canvas_p6.draw()

            def click_pG():
                plt.rcParams["font.family"] = graph_font

                SUM=data_chouten_zenbu_p["CT1"]+data_chouten_zenbu_p["CT2"]+data_chouten_zenbu_p["CT3"]
                pX = data_chouten_zenbu_p["CT1"]/SUM*100
                pY = data_chouten_zenbu_p["CT2"]/SUM*100
                pZ = data_chouten_zenbu_p["CT3"]/SUM*100
                data_concat_p = pd.concat([pX,pY,pZ,SUM],axis=1)
                data_concat_p.columns = ['X','Y','Z','SUM']

                if size_hosei_p == "":
                    data_z_SUM_p = data_choutensize_zentai
                    data_z_x_p = data_choutensize_zentai #いらない列
                elif size_hosei_p == "sizefactor":
                    if tate_iso_p == "all":
                        data_z_SUM_p=(data_choutensize_zentai*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai #いらない列
                    elif tate_iso_p == "A":
                        data_z_SUM_p=(data_choutensize_kijun*size_factor_p)**(1/3)
                        data_z_x_p=data_choutensize_zentai#いらない列
                    else:
                        messagebox.showerror("エラー","エラー")
                else:
                    messagebox.showerror('エラー','エラー')

                data_z_p = pd.concat([data_z_SUM_p,data_z_x_p],axis=1)
                data_z_p.columns = ['SUM','X']         

                zmax_p = data_z_p['SUM'].max()

                amari_p = zmax_p%10

                if amari_p == 0:
                   baisuu_p = zmax_p
                else:
                   baisuu_p = zmax_p-amari_p+10

                if textBox_p11.get() == "":
                    tB_p11 = "0"
                else:
                    tB_p11 = textBox_p11.get()

                if textBox_p12.get() == "":
                    tB_p12 = baisuu_p
                else:
                    tB_p12 = textBox_p12.get()
                global triangle_plot
                triangle_plot = 1
                global data_concat0_tri_g
                if tB_p11=="0" and tB_p12 == baisuu_p:
                   data_concat0_tri_g_p = data_concat_p
                   data_gattai0_g_p = data_chouten_zenbu_p
                else:
                    data_concat0_tri_g_p = data_concat_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]
                    data_gattai0_g_p = data_chouten_zenbu_p[(data_z_p["SUM"]>=int(tB_p11))&(data_z_p["SUM"]<=int(tB_p12))]

                fig_p=plt.figure(figsize=(7.5,7.5))
                ax_p=fig_p.add_subplot(111,projection='3d')
                ax_p.set_xticks([])
                ax_p.set_yticks([])
                ax_p.set_zticks([])
                plt.axis('off')
                ax_p.view_init(elev=12,azim=-69)

                h = np.sqrt(3.0)*0.5

                hi_list_p = []
                for i in range(0,10):
                   for j in range(0,10-i):
                      data_concat_hi_p1 = data_concat0_tri_g_p[(data_concat0_tri_g_p["X"]>=i*10)&(data_concat0_tri_g_p["X"]<=(i+1)*10)&(data_concat0_tri_g_p["Y"]>=(9-i-j)*10)&(data_concat0_tri_g_p["Y"]<=(10-i-j)*10)&(data_concat0_tri_g_p["Z"]>=j*10)&(data_concat0_tri_g_p["Z"]<=(j+1)*10)]
                      hi_p1=len(data_concat_hi_p1)
                      hi_list_p.append(hi_p1)

                for i in range(0,10):
                    for j in range(0,10-i-1):
                      data_concat_hi_p2 = data_concat0_tri_g_p[(data_concat0_tri_g_p["X"]>=i*10)&(data_concat0_tri_g_p["X"]<=(i+1)*10)&(data_concat0_tri_g_p["Y"]>=(8-i-j)*10)&(data_concat0_tri_g_p["Y"]<=(9-i-j)*10)&(data_concat0_tri_g_p["Z"]>=j*10)&(data_concat0_tri_g_p["Z"]<=(j+1)*10)]
                      hi_p2=len(data_concat_hi_p2)
                      hi_list_p.append(hi_p2)

                hi_max_p = max(hi_list_p)
                print(hi_max_p)

                for i in range(0,10):
                    for j in range(0,10-i):
                      data_concat_hi_p = data_concat0_tri_g_p[(data_concat0_tri_g_p["X"]>=i*10)&(data_concat0_tri_g_p["X"]<=(i+1)*10)&(data_concat0_tri_g_p["Y"]>=(9-i-j)*10)&(data_concat0_tri_g_p["Y"]<=(10-i-j)*10)&(data_concat0_tri_g_p["Z"]>=j*10)&(data_concat0_tri_g_p["Z"]<=(j+1)*10)]
                      hi_p=len(data_concat_hi_p)
                      if hi_p >= hi_bord_1:
                          hi_color_p = "red"
                      elif hi_p >= hi_bord_2:
                          hi_color_p = "darkorange"
                      elif hi_p >= hi_bord_3:
                          hi_color_p = "lawngreen"
                      elif hi_p >= hi_bord_4:
                          hi_color_p = "green"
                      elif hi_p >= hi_bord_5:
                          hi_color_p = "blue"
                      elif hi_p > 0:
                          hi_color_p = "lightskyblue"
                      else:
                          hi_color_p = "#ffffff00"
                      x1=[i/20+j/10,i/20+(j+1)/10,i/20+(j+1)/10,i/20+j/10]
                      y1=[i*h/10,i*h/10,i*h/10,i*h/10]
                      z1=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly1=list(zip(x1,y1,z1))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly1],facecolors=hi_color_p,linewidth=1,alpha=0.2))
                      x2=[i/20+(j+1)/10,i/20+j/10+1/20,i/20+j/10+1/20,i/20+(j+1)/10]
                      y2=[i*h/10,i*h/10+h/10,i*h/10+h/10,i*h/10]
                      z2=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly2=list(zip(x2,y2,z2))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly2],color=hi_color_p,linewidth=1,alpha=0.2))
                      x3=[i/20+j/10+1/20,i/20+j/10,i/20+j/10,i/20+j/10+1/20]
                      y3=[i*h/10+h/10,i*h/10,i*h/10,i*h/10+h/10]
                      z3=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly3=list(zip(x3,y3,z3))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly3],color=hi_color_p,linewidth=1,alpha=0.2))
                      x4=[i/20+j/10,i/20+(j+1)/10,i/20+j/10+1/20]
                      y4=[i*h/10,i*h/10,i*h/10+h/10]
                      z4=[hi_p/hi_max_p,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly4=list(zip(x4,y4,z4))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly4],color=hi_color_p,linewidth=1,alpha=0.2))
                      x5=[i/20+j/10,i/20+(j+1)/10,i/20+j/10+1/20]
                      y5=[i*h/10,i*h/10,i*h/10+h/10]
                      z5=[0,0,0]
                      poly5=list(zip(x5,y5,z5))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly5],color=hi_color_p,linewidth=1,alpha=0.2))

                for i in range(0,10):
                    for j in range(0,10-i-1):
                      data_concat_hi_p = data_concat0_tri_g_p[(data_concat0_tri_g_p["X"]>=i*10)&(data_concat0_tri_g_p["X"]<=(i+1)*10)&(data_concat0_tri_g_p["Y"]>=(8-i-j)*10)&(data_concat0_tri_g_p["Y"]<=(9-i-j)*10)&(data_concat0_tri_g_p["Z"]>=j*10)&(data_concat0_tri_g_p["Z"]<=(j+1)*10)]
                      hi_p=len(data_concat_hi_p)
                      if hi_p >= hi_bord_1:
                          hi_color_p = "red"
                      elif hi_p >= hi_bord_2:
                          hi_color_p = "darkorange"
                      elif hi_p >= hi_bord_3:
                          hi_color_p = "lawngreen"
                      elif hi_p >= hi_bord_4:
                          hi_color_p = "green"
                      elif hi_p >= hi_bord_5:
                          hi_color_p = "blue"
                      elif hi_p > 0:
                          hi_color_p = "lightskyblue"
                      else:
                          hi_color_p = "#ffffff00"
                      x1=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+(j+1)/10,i/20+j/10+3/20]
                      y1=[i*h/10+h/10,i*h/10,i*h/10,i*h/10+h/10]
                      z1=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly1=list(zip(x1,y1,z1))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly1],facecolors=hi_color_p,linewidth=1,alpha=0.2))
                      x2=[i/20+(j+1)/10,i/20+j/10+1/20,i/20+j/10+1/20,i/20+(j+1)/10]
                      y2=[i*h/10,i*h/10+h/10,i*h/10+h/10,i*h/10]
                      z2=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly2=list(zip(x2,y2,z2))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly2],color=hi_color_p,linewidth=1,alpha=0.2))
                      x3=[i/20+j/10+1/20,i/20+j/10+3/20,i/20+j/10+3/20,i/20+j/10+1/20]
                      y3=[i*h/10+h/10,i*h/10+h/10,i*h/10+h/10,i*h/10+h/10]
                      z3=[0,0,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly3=list(zip(x3,y3,z3))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly3],color=hi_color_p,linewidth=1,alpha=0.2))
                      x4=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+j/10+1/20]
                      y4=[i*h/10+h/10,i*h/10,i*h/10+h/10]
                      z4=[hi_p/hi_max_p,hi_p/hi_max_p,hi_p/hi_max_p]
                      poly4=list(zip(x4,y4,z4))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly4],color=hi_color_p,linewidth=1,alpha=0.2))
                      x5=[i/20+j/10+3/20,i/20+(j+1)/10,i/20+j/10+1/20]
                      y5=[i*h/10+h/10,i*h/10,i*h/10+h/10]
                      z5=[0,0,0]
                      poly5=list(zip(x5,y5,z5))
                      ax_p.add_collection3d(art3d.Poly3DCollection([poly5],color=hi_color_p,linewidth=1,alpha=0.2))

                for i in range(1,10):
                    ax_p.plot([2*i/40.0, 1.0-2*i/40.0],[h*2*i/20.0, h*2*i/20.0],[0,0],color='gray', lw=0.5)
                    ax_p.plot([2*i/40.0, 2*i/20.0],[h*2*i/20.0, 0.0],[0,0], color='gray', lw=0.5)
                    ax_p.plot([0.5+2*i/40.0, 2*i/20.0],[h*(1.0-2*i/20.0), 0.0],[0,0], color='gray', lw=0.5)

                ax_p.plot([0.0, 1.0],[0.0, 0.0],[0,0], color=color_axes, lw=2)
                ax_p.plot([0.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p.plot([1.0, 0.5],[0.0, h],[0,0], color=color_axes, lw=2)
                ax_p.text(0.5, h+0.15,-0.1, combobox_p11.get(), fontsize=20,ha="center", color=color_axes)
                ax_p.text(-0.15*h, -0.15/2,-0.1, combobox_p14.get(), fontsize=20,ha="center", color=color_axes)
                ax_p.text(1+0.15*h, -0.15/2,-0.1, combobox_p17.get(), fontsize=20,ha="center", color=color_axes)
                ax_p.set_xlim(0,1)
                ax_p.set_ylim(0,1)
                ax_p.set_zlim(0,1)

                fig_p.show()

            def onclick_p(event):
                print("event.button=%d,event.xdata=%f,event.ydata=%f" %(event.button,event.xdata,event.ydata))
                eventx_min = event.xdata-0.02
                eventx_max = event.xdata+0.02
                eventy_min = event.ydata-0.02
                eventy_max = event.ydata+0.02
                data_click_kouho_1=data_click_p[(data_click_p["x"]>=eventx_min)&(data_click_p["x"]<=eventx_max)]
                data_click_kouho_2=data_click_kouho_1[(data_click_kouho_1["y"]>=eventy_min)&(data_click_kouho_1["y"]<=eventy_max)]
                data_y_sa = abs(data_click_kouho_2["y"]-event.ydata)
                data_click_satuki = pd.concat([data_click_kouho_2,data_y_sa],axis=1)
                data_click_satuki.columns = ["a'","b'","c'","sum","x","y","sa"]
                data_click_minjyun = data_click_satuki.sort_values("sa")
                data_click_minjyun_r = data_click_minjyun.reset_index()
                data_click_target_gyou = data_click_minjyun_r[0:1]

                ax_p1_5.cla()

                plt.rcParams["font.family"] = graph_font

                fig_p1_5.set_facecolor(color_fig)
                ax_p1_5.set_facecolor(color_fig)

                ax_p1_5.spines['top'].set_color(color_axes)
                ax_p1_5.spines['bottom'].set_color(color_axes)
                ax_p1_5.spines['left'].set_color(color_axes)
                ax_p1_5.spines['right'].set_color(color_axes)
                ax_p1_5.tick_params(colors=color_axes)

                NP_number = int(data_click_target_gyou["index"])
                data_skip_column = int(data_large.columns.get_loc("skip"))
                element_number = int(data_skip_column)-2
                data_plot = data_large.iloc[:,int(data_skip_column):len(data_large)]
                data_time = data_large["Time"][NP_number]
                data_time_list = re.findall(r"\d+\.\d*",data_time)
                data_time_list_float = list(map(float,data_time_list))
                data_large_legend = data_large.columns[1:data_skip_column-1]

                def calc_002(n):
                    return n*0.002
                
                if iso_spe == "":
                    for i in range(0,element_number):
                        data_element_i = data_large.iloc[NP_number,data_skip_column+i+1]#[data_skip_column+i+1][NP_number]
                        data_element_list_i = re.findall(r"[+-]?[0-9]+\.[0-9]*[e]?[+-]?[0-9]*",data_element_i)
                        data_element_list_float_ini_i = list(map(float,data_element_list_i))
                        data_element_list_float_i = list(map(calc_002,data_element_list_float_ini_i))
                        ax_p1_5.plot(data_time_list_float,data_element_list_float_i,lw=2)
                elif iso_spe == "specified":
                    for i in range(0,element_number_spe):
                        data_element_i = data_large.iloc[NP_number,data_skip_column+data_large_legend_spe_index[i]+1]
                        data_element_list_i = re.findall(r"[+-]?[0-9]+\.[0-9]*[e]?[+-]?[0-9]*",data_element_i)
                        data_element_list_float_ini_i = list(map(float,data_element_list_i))
                        data_element_list_float_i = list(map(calc_002,data_element_list_float_ini_i))
                        ax_p1_5.plot(data_time_list_float,data_element_list_float_i,lw=2)


                ax_p1_5.set_xlabel("Time(s)",fontname=graph_font,color=color_axes,weight=weight_font)
                ax_p1_5.set_ylabel("Intensity(Counts)",fontname=graph_font,color=color_axes,weight=weight_font)
                if iso_spe == "":
                    ax_p1_5.legend(data_large_legend,bbox_to_anchor=(1.05,1),loc="upper left",borderaxespad=0)
                elif iso_spe == "specified":
                    ax_p1_5.legend(data_large_legend_spe,bbox_to_anchor=(1.05,1),loc="upper left",borderaxespad=0)

                canvas_p1_5.draw()

                a_dash = data_click_target_gyou["a'"]*data_click_target_gyou["sum"]/100
                b_dash = data_click_target_gyou["b'"]*data_click_target_gyou["sum"]/100
                c_dash = data_click_target_gyou["c'"]*data_click_target_gyou["sum"]/100
                a_x = a_dash*float(textBox30.get())
                b_y = b_dash*float(textBox31.get())
                c_z = c_dash*float(textBox32.get())
                a_warux = a_dash/float(textBox30.get())
                b_waruy = b_dash/float(textBox31.get())
                c_waruz = c_dash/float(textBox32.get())
                sa_per = (((1/a_x)+(a_warux+b_waruy+c_waruz)*(1/(a_dash+b_dash+c_dash))**2)**(1/2))*data_click_target_gyou["a'"]
                sb_per = (((1/b_y)+(a_warux+b_waruy+c_waruz)*(1/(a_dash+b_dash+c_dash))**2)**(1/2))*data_click_target_gyou["b'"]
                sc_per = (((1/c_z)+(a_warux+b_waruy+c_waruz)*(1/(a_dash+b_dash+c_dash))**2)**(1/2))*data_click_target_gyou["c'"]
                print("1SD_a")
                print(sa_per)
                print("1SD_b")
                print(sb_per)
                print("1SD_c")
                print(sc_per)

                ax_p3.cla()
                fig_p3.set_facecolor(color_fig)
                ax_p3.set_facecolor(color_fig)

                h=np.sqrt(3.0)*0.5

                for i in range(1,10):
                    ax_p3.plot([i/20.0, 1.0-i/20.0],[h*i/10.0, h*i/10.0], linestyle='dashed',color='gray', lw=0.5)
                    ax_p3.plot([i/20.0, i/10.0],[h*i/10.0, 0.0], linestyle='dashed',color='gray', lw=0.5)
                    ax_p3.plot([0.5+i/20.0, i/10.0],[h*(1.0-i/10.0), 0.0], linestyle='dashed',color='gray', lw=0.5)

                ax_p3.plot([0.0, 1.0],[0.0, 0.0], color=color_axes, lw=2)
                ax_p3.plot([0.0, 0.5],[0.0, h], color=color_axes, lw=2)
                ax_p3.plot([1.0, 0.5],[0.0, h], color=color_axes, lw=2)
  
                ax_p3.text(0.455, h+0.0283, combobox_p11.get(), fontsize=22, color=color_axes)
                ax_p3.text(-0.1, -0.02, combobox_p14.get(), fontsize=22, color=color_axes)
                ax_p3.text(1.02, -0.02, combobox_p17.get(), fontsize=22, color=color_axes)
                

                for i in range(1,10):
                    ax_p3.text(0.5+(10-i)/20.0+0.016, h*(1.0-(10-i)/10.0), '%d0' % i, fontsize=17, color=color_axes)
                    ax_p3.text((10-i)/20.0-0.082, h*(10-i)/10.0, '%d0' % i, fontsize=17, color=color_axes)
                    ax_p3.text(i/10.0-0.03, -0.06, '%d0' % i, fontsize=17, color=color_axes)

                ax_p3.text(-0.15,1,"Number of particles:"+str(len(data_plot_p3.dropna())),fontsize=14, color=color_axes)

                ax_p3.scatter(data_plot_p3["x3"],data_plot_p3["y3"],c=color_plot,alpha=alpha_plot,s=size_plot)

                ax_p3.scatter(data_click_target_gyou["x"],data_click_target_gyou["y"],c="orange",s=size_plot)

                target_x = data_click_target_gyou["x"]
                target_y = data_click_target_gyou["y"]

                target_y_a_1 = target_y-sa_per*h/100
                target_y_a_2 = target_y+sa_per*h/100
                target_x_a_1 = target_x+sa_per*h/100/h/2
                target_x_a_2 = target_x-sa_per*h/100/h/2

                target_x_b_1 = target_x-sb_per*h/100/h/2
                target_x_b_2 = target_x+sb_per*h/100/h/2
                target_y_b_1 = target_y-sb_per*h/100
                target_y_b_2 = target_y+sb_per*h/100

                target_x_c_1 = target_x-sc_per*h/100/h
                target_x_c_2 = target_x+sc_per*h/100/h
                target_y_c_1 = target_y
                target_y_c_2 = target_y
                
                ax_p3.plot([target_x_a_1,target_x_a_2],[target_y_a_1,target_y_a_2],color="orange",alpha=0.6,lw=1)
                ax_p3.plot([target_x_b_1,target_x_b_2],[target_y_b_1,target_y_b_2],color="orange",alpha=0.6,lw=1)
                ax_p3.plot([target_x_c_1,target_x_c_2],[target_y_c_1,target_y_c_2],color="orange",alpha=0.6,lw=1)

                canvas_p3.draw()

            fig_p3.canvas.mpl_connect("button_press_event",onclick_p)

            def savefig_p3():
                file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
                print(file_path)

                if len(file_path) != 0:
                    fig_p3.savefig(file_path)
                
            def savefig_p5():
                file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
                print(file_path)
                
                if len(file_path) != 0:
                    fig_p5.savefig(file_path)

            def savefig_p6():
                file_path = tkinter.filedialog.asksaveasfilename(defaultextension="png",filetypes=[("PNG(*.png)","*.png")])
                print(file_path)
                
                if len(file_path) != 0:
                    fig_p6.savefig(file_path)
            
            def save_csv_large_p():
                if sample_label_p == "":
                    messagebox.showerror("エラー","三角グラフを表示してください")
                else:
                    print(sample_label_p)
                    newwindow_large_p = tk.Toplevel(root)
                    newwindow_large_p.attributes("-topmost",True)
                    newwindow_large_p.geometry("500x120")
                    newwindow_large_p.title(u"範囲指定csv出力(large)")
                    label_large_p_1 = tk.Label(newwindow_large_p,text="ファイル保存",font=(label_font_jp,10),anchor=tk.W)
                    label_large_p_1.place(width=80,height=20,x=10,y=10)
                    label_large_p_2 = tk.Label(newwindow_large_p,text=sample_label_p,font=(label_font_jp,10),anchor=tk.W)
                    label_large_p_2.place(width=480,height=20,x=20,y=30)
                    label_large_p_3 = tk.Label(newwindow_large_p,text="_",font=(label_font_jp,10),anchor=tk.W)
                    label_large_p_3.place(width=10,height=20,x=20,y=50)
                    textBox_large_p_1 = ttk.Entry(newwindow_large_p)
                    textBox_large_p_1.place(width=100,height=20,x=28,y=50)
                    label_large_p_4 = tk.Label(newwindow_large_p,text="_1_NP_events_large.csv",font=(label_font_jp,10),anchor=tk.W)
                    label_large_p_4.place(width=200,height=20,x=130,y=50)

                    def OK_large_p():
                        filename_large_p = textBox_large_p_1.get()
                        data_concat_syuturyoku_hanni_large_p.to_csv("{}_{}_1_NP_events_large.csv".format(sample_label_p,filename_large_p))
                        newwindow_large_p.destroy()

                    def Cancel_large_p():
                        newwindow_large_p.destroy()

                    btn_large_p_1=tk.Button(newwindow_large_p,text="決定",command=OK_large_p)
                    btn_large_p_1.place(width=80,height=25,x=310,y=80)
                    btn_large_p_2=tk.Button(newwindow_large_p,text="キャンセル",command=Cancel_large_p)
                    btn_large_p_2.place(width=80,height=25,x=400,y=80)


            def count_size_p():
                newwindow_p7 = tk.Toplevel(newwindow_p2)
                newwindow_p7.geometry("250x200")
                newwindow_p7.title(u"カウント数 → 粒径")

                label_p68 = tk.Label(newwindow_p7,text="Size Factor = ",font=(label_font_jp,9),anchor=tk.W)
                label_p68.place(width=100,height=20,x=20,y=30)

                textBox_p64 = ttk.Entry(newwindow_p7)
                textBox_p64.place(width=50,height=20,x=120,y=30)
                textBox_p64.insert(tkinter.END,size_factor_p)    

                bln_p7_3 = tkinter.BooleanVar()
                if tate_iso_p == "A":
                    bln_p7_3.set(True)
                elif tate_iso_p == "all":
                    bln_p7_3.set(False)
                else:
                    print("checkbox_error")   
                check_p7_3 = tk.Checkbutton(newwindow_p7,variable=bln_p7_3,text="{}の大きさのみにする".format(size_kijun),font=(label_font_jp,9))
                check_p7_3.place(x=20,y=80)

                def OK_p7():
                    global size_hosei_p
                    global size_factor_p
                    global tate_iso_p
                    if textBox_p64.get() == "":
                        messagebox.showerror("エラー","Size Factorを入力してください")
                    elif float(textBox_p64.get()) <= 0:
                        messagebox.showerror("エラー","Size Factorは正の値を入力してください")
                    else:
                        size_hosei_p = "sizefactor"
                        size_factor_p = float(textBox_p64.get())
                        if bln_p7_3.get():
                            tate_iso_p = "A"
                        else:
                            tate_iso_p = "all"      
                        newwindow_p7.destroy()
            
                def Cancel_p7():
                    global size_hosei_p
                    size_hosei_p = ""
                    newwindow_p7.destroy()
            
                btn22=tk.Button(newwindow_p7,text="決定",command=OK_p7)
                btn22.place(width=80,height=25,x=20,y=140)
                btn23=tk.Button(newwindow_p7,text="キャンセル",command=Cancel_p7)
                btn23.place(width=80,height=25,x=110,y=140)
            
            menubar_p2 = tk.Menu(newwindow_p2)
            newwindow_p2.config(menu = menubar_p2)

            savefig_menu_p2 = tk.Menu(menubar_p2,tearoff = 0)
            menubar_p2.add_cascade(label = "グラフ保存",menu = savefig_menu_p2)

            tool_menu_p2 = tk.Menu(menubar_p2,tearoff = 0)
            menubar_p2.add_cascade(label = "ツール",menu = tool_menu_p2)

            savefig_menu_p2.add_command(label = "三角グラフ",command = savefig_p3)
            savefig_menu_p2.add_command(label = "組成分布",command = savefig_p5)
            savefig_menu_p2.add_command(label = "カウント分布",command = savefig_p6)
            savefig_menu_p2.add_command(label = "範囲指定csv出力(large)",command = save_csv_large_p)

            tool_menu_p2.add_command(label = "カウント数 → 粒径",command = count_size_p)

            btn_p3 = tk.Button(newwindow_p2,text="読み込み",command=click_p2,font=(label_font_jp,8))
            btn_p3.place(width=160,height=25,x=100,y=200)
            btn_p4 = tk.Button(newwindow_p2,text="表示",command=click_pC,font=(label_font_jp,8))
            btn_p4.place(width=50,height=25,x=215,y=359)
            btn_p5 = tk.Button(newwindow_p2,text="全体表示",command=click_pA,font=(label_font_jp,9))
            btn_p5.place(width=90,height=25,x=60,y=440)
            btn_p6=tk.Button(newwindow_p2,text="表示",command=click_pD,font=(label_font_jp,9))
            btn_p6.place(width=50,height=25,x=215,y=469)
            btn_p7=tk.Button(newwindow_p2,text="表示",command=click_pE,font=(label_font_jp,9))
            btn_p7.place(width=50,height=25,x=215,y=609)
            btn_p8=tk.Button(newwindow_p2,text="表示",command=click_pF,font=(label_font_jp,9))
            btn_p8.place(width=50,height=25,x=215,y=809)
            btn_p9=tk.Button(newwindow_p2,text="表示",command=click_pG,font=(label_font_jp,9))
            btn_p9.place(width=50,height=25,x=180,y=509)

            textBox_p11 = tk.Entry(newwindow_p2)
            textBox_p11.place(width=40,height=20,x=110,y=360)
            textBox_p12 = tk.Entry(newwindow_p2)
            textBox_p12.place(width=40,height=20,x=160,y=360)
            label_p23=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font,10))
            label_p23.place(width=10,height=20,x=150,y=360)
            label_p24=tk.Label(newwindow_p2,text="カウント",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p24.place(width=50,height=20,x=60,y=360)       

            textBox_p13=tk.Entry(newwindow_p2)
            textBox_p13.place(width=40,height=20,x=110,y=470)
            textBox_p14=tk.Entry(newwindow_p2)
            textBox_p14.place(width=40,height=20,x=160,y=470)
            label_p25=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font,10))
            label_p25.place(width=10,height=20,x=150,y=470)
            label_p26=tk.Label(newwindow_p2,text="カウント",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p26.place(width=50,height=20,x=60,y=470)      

            textBox_p15=tk.Entry(newwindow_p2)
            textBox_p15.place(width=40,height=20,x=156,y=550)
            textBox_p16=tk.Entry(newwindow_p2)
            textBox_p16.place(width=40,height=20,x=205,y=550)
            combobox_p20=ttk.Combobox(newwindow_p2,state='readonly',values=[1,2,4,5,10,20,25,50,100])
            combobox_p20.current(4)      
            combobox_p20.place(width=40,height=20,x=110,y=610)
            label_p27=tk.Label(newwindow_p2,text=":",bg=color8,font=(label_font_jp,10))
            label_p27.place(width=3,height=20,x=151,y=550)
            label_p28=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font_jp,10))
            label_p28.place(width=10,height=20,x=195,y=550)
            label_p29=tk.Label(newwindow_p2,text="%",bg=color8,font=(label_font_jp,9))
            label_p29.place(width=20,height=20,x=248,y=550) 
            label_p30=tk.Label(newwindow_p2,text="固定",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p30.place(width=50,height=20,x=60,y=550) 
            label_p31=tk.Label(newwindow_p2,text="表示",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p31.place(width=50,height=20,x=60,y=580) 
            label_p32=tk.Label(newwindow_p2,text="幅",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p32.place(width=50,height=20,x=60,y=610) 

            textBox_p17=tk.Entry(newwindow_p2)
            textBox_p17.place(width=40,height=20,x=110,y=690)
            textBox_p18=tk.Entry(newwindow_p2)
            textBox_p18.place(width=40,height=20,x=160,y=690)
            textBox_p19=tk.Entry(newwindow_p2)
            textBox_p19.place(width=40,height=20,x=110,y=720)
            textBox_p20=tk.Entry(newwindow_p2)
            textBox_p20.place(width=40,height=20,x=160,y=720)
            textBox_p21=tk.Entry(newwindow_p2)
            textBox_p21.place(width=40,height=20,x=110,y=750)
            textBox_p22=tk.Entry(newwindow_p2)
            textBox_p22.place(width=40,height=20,x=160,y=750)
            label_p33=tk.Label(newwindow_p2,text="A"+":",bg=color8,font=(label_font_jp,10))
            label_p33.place(width=20,height=20,x=80,y=690)
            label_p34=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font_jp,10))
            label_p34.place(width=10,height=20,x=150,y=690)
            label_p35=tk.Label(newwindow_p2,text="%",bg=color8,font=(label_font_jp,9))
            label_p35.place(width=20,height=20,x=203,y=690)
            label_p36=tk.Label(newwindow_p2,text="B"+":",bg=color8,font=(label_font_jp,10))
            label_p36.place(width=20,height=20,x=80,y=720)
            label_p37=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font_jp,10))
            label_p37.place(width=10,height=20,x=150,y=720)
            label_p38=tk.Label(newwindow_p2,text="%",bg=color8,font=(label_font_jp,9))
            label_p38.place(width=20,height=20,x=203,y=720)
            label_p39=tk.Label(newwindow_p2,text="C"+":",bg=color8,font=(label_font_jp,10))
            label_p39.place(width=20,height=20,x=80,y=750)
            label_p40=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font_jp,10))
            label_p40.place(width=10,height=20,x=150,y=750)
            label_p41=tk.Label(newwindow_p2,text="%",bg=color8,font=(label_font_jp,9))
            label_p41.place(width=20,height=20,x=203,y=750)
            textBox_p23=tk.Entry(newwindow_p2)
            textBox_p23.place(width=40,height=20,x=110,y=780)
            textBox_p24=tk.Entry(newwindow_p2)
            textBox_p24.place(width=40,height=20,x=160,y=780)      
            label_p42=tk.Label(newwindow_p2,text="-",bg=color8,font=(label_font,10))
            label_p42.place(width=10,height=20,x=150,y=780)
            label_p43=tk.Label(newwindow_p2,text="カウント",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p43.place(width=50,height=20,x=60,y=780) 
            textBox_p25=tk.Entry(newwindow_p2)
            textBox_p25.place(width=40,height=20,x=110,y=810) 
            label_p44=tk.Label(newwindow_p2,text="幅",anchor=tk.CENTER,bg=color8,font=(label_font_jp,9))
            label_p44.place(width=50,height=20,x=60,y=810) 

    def Cancelp():
        newwindow_p.destroy()

    btn_p1=tk.Button(newwindow_p,text="決定",command=OKp)
    btn_p1.place(width=80,height=25,x=680,y=230)
    btn_p2=tk.Button(newwindow_p,text="キャンセル",command=Cancelp)
    btn_p2.place(width=80,height=25,x=780,y=230)

menubar = tk.Menu(root)
root.config(menu = menubar)

save_menu = tk.Menu(menubar, tearoff = 0)
menubar.add_cascade(label = 'ファイル', menu = save_menu)

savefig_menu = tk.Menu(menubar, tearoff = 0)
menubar.add_cascade(label = 'グラフ保存', menu = savefig_menu)

tool_menu = tk.Menu(menubar, tearoff = 0)
menubar.add_cascade(label = 'ツール', menu = tool_menu)

setting_menu = tk.Menu(menubar, tearoff = 0)
menubar.add_cascade(label = '設定', menu = setting_menu)

savefig_menu.add_command(label = '三角グラフ',command = savefig3)
savefig_menu.add_command(label = '組成分布',command = savefig5)
savefig_menu.add_command(label = 'カウント分布',command = savefig6)
savefig_menu.add_command(label = '範囲指定csv出力(large)',command = save_csv_large)
save_menu.add_command(label = '名前を付けて保存',command = savelist)
save_menu.add_command(label = '開く',command = openlist)
save_menu.add_separator()
save_menu.add_command(label = '終了',command = end)
setting_menu.add_command(label = 'レイアウト設定',command = layout)
setting_menu.add_command(label = '三次元組成分布設定',command = distribution3d)
tool_menu.add_command(label = '表示同位体指定',command = isospe)
tool_menu.add_command(label = '参照値',command = sansyou)
tool_menu.add_command(label = "頂点設定",command = choutenplus)
tool_menu.add_command(label = "カウント数 → 粒径",command = count_size)

root.mainloop()