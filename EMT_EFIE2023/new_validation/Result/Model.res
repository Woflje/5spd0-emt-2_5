MWS Result File Version 20150206
size=i:30

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:E-field,Theta=0.rd1

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:RefSpectrum_pw.sig

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:e-field (f=15)_pw.m3d

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:h-field (f=15)_pw.m3d

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:solverstart
result=s:0
files=s:PBAConnectivity.axg

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:solverstart
result=s:0
files=s:PBAMeshFeedback.axg

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:World.fid

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:survivemeshadapt
result=s:1
files=s:model.gex

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:survivemeshadapt
result=s:1
files=s:PP.sid

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:survivemeshadapt
result=s:1
files=s:PP.fmm

type=s:HIDDENITEM
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
files=s:ml_info.dat

type=s:FOLDER
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:persistent
result=s:0
treepath=s:1D Results

type=s:MESH_FEEDBACK
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:solverstart
result=s:0
treepath=s:Mesh\Information\PBA
files=s:PBAMeshFeedback.rex
ylabel=s:Mesh Feedback

type=s:MESH_FEEDBACK
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:solverstart
result=s:0
treepath=s:Mesh\Information\Connectivity
files=s:PBAConnectivity.rex
ylabel=s:Mesh Feedback

type=s:XYSIGNAL2
subtype=s:user
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\Port signals\Plane wave
files=s:plw.sig

type=s:HFIELD3D
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:2D/3D Results\H-Field\h-field (f=15) [pw]
files=s:h-field (f=15)_pw.m3d
files=s:h-field (f=15)_pw_m3d.rex

type=s:SURFACECURRENT
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:2D/3D Results\Surface Current\surface current (f=15) [pw]
files=s:h-field (f=15)_pw.m3d
files=s:h-field (f=15)_pw_m3d_sct.rex

type=s:EFIELD3D
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:2D/3D Results\E-Field\e-field (f=15) [pw]
files=s:e-field (f=15)_pw.m3d
files=s:e-field (f=15)_pw_m3d.rex

type=s:XYSIGNAL2
subtype=s:energy
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\Energy\Energy [pw]
files=s:pw.eng

type=s:FARFIELD
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:Farfields\farfield (f=15) [pw]
files=s:farfield (f=15)_pw.ffm
ylabel=s:farfield (f=15) [pw]

type=s:XYSIGNAL2
subtype=s:complex
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\Power\Excitation [pw]\Power Scattered
files=s:FarfieldMetaData_pw_RadPower.sig

type=s:XYSIGNAL2
subtype=s:linear
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\Cross Sections\Total RCS [pw]
files=s:FarfieldMetaData_pw_TotRCS.sig

type=s:XYSIGNAL2
subtype=s:linear
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\Cross Sections\Total ACS [pw]
files=s:FarfieldMetaData_pw_TotACS.sig

type=s:FARFIELD1DCUT
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:Farfields\Farfield Cuts\Excitation [pw]\Phi=0\farfield (f=15)
files=s:Farfield_Cut_farfield (f=15)_Phi=0_[pw]_0.sig

type=s:FARFIELD1DCUT
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:Farfields\Farfield Cuts\Excitation [pw]\Phi=90\farfield (f=15)
files=s:Farfield_Cut_farfield (f=15)_Phi=90_[pw]_0.sig

type=s:FARFIELD1DCUT
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:Farfields\Farfield Cuts\Excitation [pw]\Theta=90\farfield (f=15)
files=s:Farfield_Cut_farfield (f=15)_Theta=90_[pw]_0.sig

type=s:RESULT_0D
problemclass=s::8:1000
visibility=s:hidden
creation=s:internal
lifetime=s:rebuild
result=s:1
treepath=s:1D Results\AutomaticRunInformation
files=s:AutomaticRunInformation

type=s:TABLE
subtype=s:farfield polar linear
problemclass=s::8:1000
visibility=s:visible
creation=s:internal
lifetime=s:surviveparchange
result=s:1
treepath=s:Tables\1D Results\E-field,Theta=0
files=s:E-field,Theta=0.rt1
files=s:E-field,Theta=0.rd1

type=s:XYSIGNAL2
subtype=s:user
problemclass=s::4:3
visibility=s:visible
creation=s:internal
lifetime=s:persistent
result=s:0
treepath=s:Excitation Signals\default
files=s:signal_default_lf.sig

type=s:XYSIGNAL2
subtype=s:user
problemclass=s::0:0
visibility=s:visible
creation=s:internal
lifetime=s:persistent
result=s:0
treepath=s:Excitation Signals\default
files=s:signal_default.sig

