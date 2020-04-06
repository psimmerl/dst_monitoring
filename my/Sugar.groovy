package my

import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase

class Sugar {
  static void enable() {
    H2F.metaClass.projectionX = {name,title=""->
      H1F h1 = delegate.projectionX()
      h1.setName(name)
      h1.setTitle(title)
      return h1
    }

    H2F.metaClass.projectionY = {name,title=""->
      H1F h1 = delegate.projectionY()
      h1.setName(name)
      h1.setTitle(title)
      return h1
    }

    TDirectory.metaClass.writeDataSet = {
      def pwd = delegate.getDir().getDirectoryPath()
      it.getName().split('/').dropRight(1).findAll().each{
        delegate.mkdir(it)
        delegate.cd(it)
      }
      delegate.add(it.getName().split('/')[-1], it)
      if(pwd=='/') delegate.cd()
      else delegate.cd(pwd)
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

