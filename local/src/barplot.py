import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 1+5:
	sys.exit("""
	Usage: %s data_file title label1 label2 output_image [output_image2 ...]

	'data_file': sample file1 file2 clones1 clones2
""" % sys.argv[0])

DATA = sys.argv[1]
TITLE = sys.argv[2]
LABEL1 = sys.argv[3]
LABEL2 = sys.argv[4]
OUT = sys.argv[5:]

data = [i.strip('\n').split('\t') for i in open(DATA)]

xlabels = [i[0] for i in data]

X = range(len(data))
Y1 = [int(i[3]) for i in data]
Y2 = [-int(i[4]) for i in data]

YTICKS = range(min(Y2),max(Y1)+1)

fig, ax = plt.subplots()
ax.bar(X, Y1, label=LABEL1, zorder=2)
ax.bar(X, Y2, label=LABEL2, zorder=2)
ax.set_xticks(X)
ax.set_xticklabels(xlabels, rotation=90)
ax.set_yticks(YTICKS)
ax.set_ylabel('N. of clones')
#plt.bar(X, Y1, label=LABEL1, zorder=2)
#plt.bar(X, Y2, label=LABEL2, zorder=2)
#plt.gca().set_yticks(YTICKS)
#plt.gca().set_xticks(X)
#plt.gca().set_xticklabels(xlabels, rotation=90)
#plt.set_ylabel('N. of clones')
plt.grid(axis = 'y')

plt.legend()
plt.title(TITLE)
plt.tight_layout()

for out in OUT:
	print("Saving %s ..." % out)
	plt.savefig(out)
