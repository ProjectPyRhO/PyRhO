# Hay2011.py

#['import3d.hoc', 'models/L5PCbiophys3.hoc', 'models/L5PCtemplate.hoc']
self.h.load_file('import3d.hoc')
self.h.load_file('models/L5PCbiophys3.hoc')
self.h.load_file('models/L5PCtemplate.hoc')
self.h('objref L5PC') #h('objref '+cell)
self.h.L5PC = self.h.L5PCtemplate("morphologies/cell1.asc")
return self.h.L5PC