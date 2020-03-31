package clas12.groovy

import org.jlab.groot.data.TDirectory
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase

class Sugar {
  static void enable() {
    TDirectory.metaClass.writeDataSet = {
      def pwd = delegate.getDir().getDirectoryPath()
      it.getName().split('/').dropRight(1).findAll().each{
        delegate.mkdir(it)
        delegate.cd(it)
      }
      delegate.add(it.getName().split('/')[-1], it)
      delegate.cd(pwd)
    }

    LorentzVector.metaClass.static.withPID = {pid,x,y,z->
      double m = PDGDatabase.getParticleMass(pid)
      double e = Math.sqrt(x*x+y*y+z*z+m*m)
      return new LorentzVector(x,y,z,e)
    }

    LorentzVector.metaClass.plus = {
      LorentzVector a = new LorentzVector(delegate)
      a.add(it)
      return a
    }

    LorentzVector.metaClass.minus = {
      LorentzVector a = new LorentzVector(delegate)
      a.sub(it)
      return a
    }

  }
}

