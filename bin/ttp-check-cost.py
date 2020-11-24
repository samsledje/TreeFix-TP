import sys
import re
import subprocess as sp
from treefix_tp.models.fitchmodel import FitchModel

try:
    oldTreFile, newTreFile = sys.argv[1:]
except ValueError:
    print('usage: python ttp-check-cost.py [old tree file] [new tree file]')
    sys.exit(0)

FM = FitchModel()

def get_cost(f):
    patt = r"Score: (?P<cost>\d+)"
    args = [FM.cmd, f]
    try:
        proc = sp.Popen(args, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True)
        ret = proc.wait()
    except:
        raise Exception("fitch.linux failed")
    if ret != 0:
        raise Exception("fitch.linux failed with returncode {}".format(ret))

    cost = None
    for line in proc.stdout:
        m = re.match(patt, line)
        if m:
            cost = int(m.group("cost"))
        break
    assert cost is not None

    return cost

old_cost = get_cost(oldTreFile)
new_cost = get_cost(newTreFile)

print("{} Cost: {}".format(oldTreFile, old_cost))
print("{} Cost: {}".format(newTreFile, new_cost))
print("Cost Decrease: {}".format(old_cost - new_cost))
