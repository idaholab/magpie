import sys

class DecayMode:
  def __init__(self, A, Z, T, E, I):
    self._A = A
    self._Z = Z
    self._T = T
    self._E = E
    self._I = I

lines = [line.rstrip('\n') for line in open(sys.argv[1])]

counter = 1
alpha_decay_map = {}
for line in lines:
  if counter != 1:
    try:
      # column 1-4 is A
      A = int(line[:5])
      # column 12-16
      Z = int(line[14:19])
      # half-life 73-92
      T = float(line[82:102])
      # alpha E value 137-149
      E = float(line[137:149])
      # intensity 185-197
      I = float(line[185:197])
      ZAID = Z * 10000 + A * 10
      if ZAID in alpha_decay_map.keys():
        alpha_decay_map[ZAID].append(DecayMode(A, Z, T, E, I))
      else:
        alpha_decay_map[ZAID] = [DecayMode(A, Z, T, E, I)]
    except:
      print "Warning. Line("+str(counter)+"):"
      print line
      print "could not be parsed"
  counter += 1

# normalize
for ZAID in alpha_decay_map.keys():
  I = 0.0
  for j in range(len(alpha_decay_map[ZAID])):
    I += alpha_decay_map[ZAID][j]._I

  for j in range(len(alpha_decay_map[ZAID])):
    alpha_decay_map[ZAID][j]._I /= I
    # we want the alpha energy in eV but it's given in keV
    alpha_decay_map[ZAID][j]._E *= 1.0e3

# print the table for reading into C++
string = "%d\n" % (len(alpha_decay_map))
for ZAID in alpha_decay_map.keys():
  string += "%d %d\n" % (ZAID, len(alpha_decay_map[ZAID]))
  for j in range(len(alpha_decay_map[ZAID])):
    string += "%d %d %25.16e %25.16e %25.16e\n" % (alpha_decay_map[ZAID][j]._Z, alpha_decay_map[ZAID][j]._A,\
                                                   alpha_decay_map[ZAID][j]._T, alpha_decay_map[ZAID][j]._E,\
                                                   alpha_decay_map[ZAID][j]._I)

print string
