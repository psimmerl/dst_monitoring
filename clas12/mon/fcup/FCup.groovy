package clas12.mon.fcup

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class FCup {
  def hists = new ConcurrentHashMap()
  def fcentry = new ConcurrentHashMap()


  def banknames = ['RUN::config', 'RUN::scaler']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (cnfb,scaler) = banknames.collect{event.getBank(it)}
      def ts = cnfb.getLong("timestamp", 0)
      fcentry[ts] = [scaler.getFloat('fcup',0), scaler.getFloat('fcupgated',0)]
    }
  }


  def finish() {
    def fc0 = null, fc1 = null, ts0 = null

    def timeline = fcentry.sort().collect{ts,fcs->[ts,*fcs]}.unique()

    def data = (0..timeline.size-2).collect{i->
      def (dt,dfc0,dfc1) = [0,1,2].collect{timeline[i+1][it]-timeline[i][it]}
      dt *= 4*1e-9
      return [dt, dfc0, dfc1, dfc0/dt, dfc1/dt]
    }

    def h0fcdt = new H1F("h0dt", "time between FC readings;time [sec]", 120,0,data.max{it[0]}[0])
    def h0fcdq0 = new H1F("h0ugtdQ", "ungated charge between FC readings:charge [nC]", 120,0,data.max{it[1]}[1])
    def h0fcdq1 = new H1F("h0gtdQ", "gated charge between FC readings:charge [nC]", 120,0,data.max{it[2]}[2])
    def h0curr0 = new H1F("h0ugtdI", "ungated current;current [nA]", 120,0,data.max{it[3]}[3])
    def h0curr1 = new H1F("h0gtdI", "gated current;current [nA]", 120,0,data.max{it[4]}[4])
    def h1fcdt = new H1F("h1dt", "time between FC readings;time [sec]", 120,0,0.5)
    def h1fcdq0 = new H1F("h1ugtdQ", "charge between FC readings:charge [nC]", 120,0,3)
    def h1fcdq1 = new H1F("h1gtdQ", "charge between FC readings:charge [nC]", 120,0,3)
    def h1curr0 = new H1F("h1ugtdI", "current;current [nA]", 120,0,70)
    def h1curr1 = new H1F("h1gtdI", "current;current [nA]", 120,0,70)

    data.each{dt,dfc0,dfc1,curr0,curr1->
      h0fcdt.fill(dt)
      h0fcdq0.fill(dfc0)
      h0fcdq1.fill(dfc1)
      h0curr0.fill(curr0)
      h0curr1.fill(curr1)

      h1fcdt.fill(dt)
      h1fcdq0.fill(dfc0)
      h1fcdq1.fill(dfc1)
      h1curr0.fill(curr0)
      h1curr1.fill(curr1)
    }

    [h0fcdt,h0fcdq0,h0fcdq1,h0curr0,h0curr1,h1fcdt,h1fcdq0,h1fcdq1,h1curr0,h1curr1].each{hists[it.getName()] = it}
  }
}
