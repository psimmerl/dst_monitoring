package clas12.mon.electron

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class ECounts {
  def hists = new ConcurrentHashMap()
  def fcentry = new ConcurrentHashMap()
  def elentry = new ConcurrentHashMap()

  def processEvent(event) {
    if(event.hasBank("RUN::config")) {
      def cnfb = event.getBank("RUN::config")
      def ts = cnfb.getLong("timestamp", 0)

      if(event.hasBank("RUN::scaler")) {
        def scaler = event.getBank("RUN::scaler")
        def fcg = scaler.getFloat('fcupgated',0)
        fcentry[ts] = fcg
      }

     if(event.hasBank("REC::Particle")) {
        def partb = event.getBank("REC::Particle")
        def iele = (0..<partb.rows()).find{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0}
        if(iele!=null) elentry[ts] = 0
      }
    }
  }


  def finish() {
   def tline = fcentry.collect{ts,fcg->[ts,0,fcg]}
   tline.addAll(elentry.collect{[it.key,1,0]})
   tline.sort{a,b->a[0]<=>b[0]?:a[1]<=>b[1]}

   def data = []
   tline.each{ts,id,fcg->
     if(id==0) data.add([ts, fcg, 0])
     else if(data) data[-1][1]++
   }
   data.each{println it}
   
/*
    def fc0 = null, fc1 = null, ts0 = null

    def timeline = fcentry.sort().collect{ts,fcs->[ts,*fcs]}

    def data = (0..timeline.size-2).collect{i->
      def (dt,dfc0,dfc1) = [0,1,2].collect{timeline[i+1][it]-timeline[i][it]}
      dt *= 4*1e-9
      return [dt, dfc0, dfc1, dfc0/dt, dfc1/dt]
    }

    def h0fcdt = new H1F("hdt_0", "time between FC readings;time [sec]", 200,0,data.max{it[0]}[0])
    def h0fcdq0 = new H1F("hugtdQ_0", "ungated charge between FC readings:charge [nC]", 300,0,data.max{it[1]}[1])
    def h0fcdq1 = new H1F("hgtdQ_0", "gated charge between FC readings:charge [nC]", 300,0,data.max{it[2]}[2])
    def h0curr0 = new H1F("hugtdI_0", "ungated current;current [nA]", 200,0,data.max{it[3]}[3])
    def h0curr1 = new H1F("hgtdI_0", "gated current;current [nA]", 200,0,data.max{it[4]}[4])
    def h1fcdt = new H1F("hdt_1", "time between FC readings;time [sec]", 200,0,0.1)
    def h1fcdq0 = new H1F("hugtdQ_1", "charge between FC readings:charge [nC]", 300,0,3)
    def h1fcdq1 = new H1F("hgtdQ_1", "charge between FC readings:charge [nC]", 300,0,3)
    def h1curr0 = new H1F("hugtdI_1", "current;current [nA]", 300,0,60)
    def h1curr1 = new H1F("hgtdI_1", "current;current [nA]", 300,0,60)

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
*/
  }
}
