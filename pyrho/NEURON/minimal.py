# minNEURON.py

cell = h.SectionList()
soma = h.Section(name='soma') #create soma
soma.push()
#h.topology()

# Geometry
soma.nseg = 1
soma.L = 20
soma.diam = 20

# Biophysics
for sec in h.allsec():
    sec.Ra = 100
    sec.cm = 1
    sec.insert('pas')
    #sec.insert('hh') #    insert hh
    cell.append(sec)

#h('objref rho')
#h('rho = new ChR(0.5)')


#h.rho.Er = Prot.phis[0]
#setattr(h.rho, 'del', Prot.pulses[0][0]) # rho.del will not work because del is reserved word in python
#h.rho.ton = Prot.onDs[0]
#h.rho.toff = Prot.offDs[0]
#h.rho.num = Prot.nPulses
#h.rho.gbar = RhO.g/20000

# Pick a rhodopsin to record from
#rhoRec = h.ChR_apic.o(7)

h.pop_section()



return cell