import os, sys
for i in range(100):
	step = i * 10
	print step
	os.system("python MyMD.py --if 1FW4_med.rvc --dt " + str(float(step)))
	os.system("python MyMD.py --if 1FW4_noCa_med.rvc --dt " + str(float(step)))
	os.system("python MyMD.py --if 1FW4_MUT_med.rvc --dt " + str(float(step)))
	print "\n\n"
