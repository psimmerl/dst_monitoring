package clas12.mon.fcup

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class FCup_mon {
  def hists = new ConcurrentHashMap()
  def fcentry = new ConcurrentHashMap()


  def banknames = ['RUN::config', 'RUN::scaler']
  def processEvent(event) {
    if(banknames.every{event.hasBank(it)}) {
      def (cnfb,scaler) = banknames.collect{event.getBank(it)}
      def ts = cnfb.getLong("timestamp", 0)
      def t0 = cnfb.getInt("unixtime", 0)
      fcentry.computeIfAbsent(t0, {[]}).add([ts, scaler.getFloat('fcup',0), scaler.getFloat('fcupgated',0)])
    }
  }


  def finish() {
    def fc0 = null, fc1 = null, ts0 = null

    def timeline = fcentry.sort().collectMany{t0,fcs->fcs.sort{a,b->a[0]<=>b[0]}.collect{[t0,*it]}}

def tl = timeline
    def data = (0..timeline.size-2).collect{i->
      def (dt,dts,dfc0,dfc1) = [0,1,2,3].collect{timeline[i+1][it]-timeline[i][it]}
      dts *= 4*1e-9
if(dt>1 || dts>1){
println("+=======================================")
println([tl[i+1][1],tl[i][1],dts,dt])
println([tl[i+1][2],tl[i][2],dfc0])
println([tl[i+1][3],tl[i][3],dfc1])
}
      return [dts, dfc0, dfc1, dfc0/dts, dfc1/dts]
    }

    def minuQ = Math.min(0,data.collect{it[1]}.min())
    def maxuQ = data.collect{it[1]}.max()
    def duQ = 0.1*(maxuQ-minuQ)
    def mingQ = Math.min(0,data.collect{it[2]}.min())
    def maxgQ = data.collect{it[2]}.max()
    def dgQ = 0.1*(maxgQ-mingQ)

    def minuI = Math.min(0,data.collect{it[3]}.min())
    def maxuI = data.collect{it[3]}.max()
    def duI = 0.1*(maxuI-minuI)
    def mingI = Math.min(0,data.collect{it[4]}.min())
    def maxgI = data.collect{it[4]}.max()
    def dgI = 0.1*(maxgI-mingI)

    def h0fcdt = new H1F("hdt_0", "time between FC readings [full range];time [sec]", 200,0,data.max{it[0]}[0])
    def h0fcdq0 = new H1F("hugtdQ_0", "ungated charge between FC readings [full range];charge [nC]", 300,minuQ-duQ,maxuQ+duQ)
    def h0fcdq1 = new H1F("hgtdQ_0", "gated charge between FC readings [full range];charge [nC]", 300,mingQ-dgQ,maxgQ+dgQ)
    def h0curr0 = new H1F("hugtdI_0", "ungated current [full range];current [nA]", 200,minuI-duI,maxuI+duI)
    def h0curr1 = new H1F("hgtdI_0", "gated current [full range];current [nA]", 200,mingI-dgI,maxgI+dgI)
    def h1fcdt = new H1F("hdt_1", "time between FC readings [fixed range];time [sec]", 200,0,0.1)
    def h1fcdq0 = new H1F("hugtdQ_1", "ungated charge between FC readings [fixed range];charge [nC]", 300,0,3)
    def h1fcdq1 = new H1F("hgtdQ_1", "gated charge between FC readings [fixed range];charge [nC]", 300,0,3)
    def h1curr0 = new H1F("hugtdI_1", "ungated current [fixed range];current [nA]", 300,0,170)
    def h1curr1 = new H1F("hgtdI_1", "gated current [fixed range];current [nA]", 300,0,170)

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
