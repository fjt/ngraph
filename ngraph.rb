#!/usr/bin/ruby
class Numeric
  def to_pi
    if self > 0
      (i = self.to_i) + (if rand < (self - i); 1; else; 0;end)
    else
      (i = self.to_i) - (if rand < (i - self); 1; else; 0;end)
    end
  end
end

def shutup(&block)
  o = STDOUT.dup; e = STDERR.dup
  STDOUT.reopen("/dev/null"); STDERR.reopen("/dev/null")
  block.call
  STDOUT.reopen(o); STDERR.reopen(e)
end

require "Nbody"

# class Nbody < Hash ## data recovering cheat code here
#   def method_missing(mtd, args=nil) 
#     if not args
#       self[mtd.to_s]
#     else
#       self[mtd.to_s.gsub("=", "")]=args
#     end
#   end
# end
 
class Ngraph
  def Ngraph.read_net(src)
    src=src.split("\n")
    vnum = src.first.strip.split(' ').last.to_i
    vertice = src[1..vnum].map{|str|str.split(' ').last.gsub('"', '')}
    linkindice = src[vnum+2..-1].map{|e|e=e.split(' ');e[0..1].map{|s|s.to_i - 1} if e.last.to_f >= 0.5}.compact
    g=Ngraph.new
    g.vertex = vertice
    g.edge = linkindice.map{|e|e.map{|i|vertice[i]}}
    g
  end

  def Ngraph.fullerene
    ret=Ngraph.new
    ret.vertex=(0..59).to_a
    ret.edge=[[2, 3], [1, 2], [0, 1], [0, 4], [2, 7], [1, 6], [0, 5], [8, 3], [3, 4], [9, 4], [22, 10], [10, 5], [19, 5], [18, 6], [17, 6], [16, 7], [15, 7], [14, 8], [13, 8], [12, 9], [11, 9], [19, 18], [17, 16], [14, 15], [13, 12], [11, 10], [20, 19], [38, 18], [36, 17], [34, 16], [32, 15], [30, 14], [28, 13], [26, 12], [24, 11], [20, 22], [38, 36], [34, 32], [30, 28], [26, 24], [21, 20], [39, 38], [37, 36], [35, 34], [33, 32], [28, 29], [30, 31], [26, 27], [24, 25], [23, 22], [23, 25], [21, 39], [37, 35], [33, 31], [29, 27], [45, 31], [46, 33], [47, 35], [48, 37], [49, 39], [40, 21], [41, 23], [42, 25], [43, 27], [44, 29], [40, 41], [49, 48], [47, 46], [45, 44], [43, 42], [51, 40], [51, 49], [59, 48], [59, 47], [57, 46], [57, 45], [55, 44], [55, 43], [53, 42], [53, 41], [52, 53], [50, 51], [58, 59], [56, 57], [54, 55], [50, 52], [50, 58], [58, 56], [56, 54], [54, 52]]
    ret
  end

  def Ngraph.platonic(faces)
    case faces
    when 4
      Ngraph.complete(4)
    when 6
      Ngraph.cube(3)
    when 8
      g=Ngraph.new
      g.vertex=(0..5).to_a
      g.edge=(0..2).map{|i|[i, (i+1)%3]} + (3..5).map{|i|[i, (i+1)%3+3]} +(0..2).map{|i|[i, (i+3)]} + (3..5).map{|i|[i, (i+1)%3]}
      g
    when 12
      g=Ngraph.new
      g.vertex=(0..19).to_a
      g.edge=(0..4).map{|i|[i, (i+1)%5]}+
        (0..4).map{|i|[i, i+5]}+
        (15..19).map{|i|[i, (i+1)%5+15]}+
        (15..19).map{|i|[i, (i-5)]}+
        [(5..9).to_a, (10..14).to_a].transpose+
        [(0..4).map{|v|(v+1)%5+5}, (10..14).to_a].transpose
      g
    when 20
      g=Ngraph.new
      g.vertex=(0..11).to_a
      g.edge=(0..4).map{|i|[0, i+1]}+
        (0..4).map{|i|[i+1, (i+1)%5+1]}+
        (0..4).map{|i|[11, i+6]}+
        (0..4).map{|i|[i+6, (i+1)%5+6]}+
        [(1..5).to_a, (6..10).to_a].transpose+
        [(0..4).map{|v|(v+4)%5+1}, (6..10).to_a].transpose
      g
    else
      nil
    end
  end

  def Ngraph.ktree(depth = 4, branch = 3)
    g = Ngraph.new
    n = (branch**depth - 1)/(branch - 1)
    a = Array.new
    n.times{|i| a.push(i)}
    g.vertex = a
    a = Array.new
    (n-1).times{|i| a.push([i/branch, i+1])}
    g.edge = a
    g
  end


  def Ngraph.cube(dim)
    cube = Ngraph.new
    cube.vertex = Array.new(2**dim){|i|i}
    nrms = Array.new(dim){|i|2**i}
    e = [];cube.vertex.each{|v|nrms.each{|f| e.push([v, (v^f)].sort)}}
    cube.edge = e.uniq
    cube
  end

  def Ngraph.erdos_reny(size:10000, prob: nil)
    lim=size-1
    prob = log(size)/size.to_f if not prob
    links=[]
    (0..lim).each{|i|
      ((i+1)..lim).each{|ii| links.push([i, ii]) if rand < prob}}
    er=Ngraph.new; er.vertex=(0..lim).to_a; er.diredge=links
    er
  end

  def Ngraph.random_graph(n, p = nil)
    ## obsolete erdos-reny model do not use
    g = Ngraph.new
    g.vertex = Array.new(n){|i|i}
    vn = g.vertex.length-1
    p = 2/g.vertex.length.to_f if p.nil?

    g.edge = g.vertex.combination(2).inject([]){|ret, e|
      if rand(0) < p
	ret+= [e]
      else
	ret
      end}
    g
  end

  def Ngraph.ba(vn, hd)
    # BA model algorithm
    hd < vn or raise "hub degree #{hd} must be less than num of vertice #{vn}."
    g=Ngraph.new
    cans = Array.new(hd){|i| i} + Array.new(hd, hd)
    eg = Array.new(hd){|i| [i, hd]}
    src = hd + 1
    while src < vn
      fg=Array.new
      hd.times{
	dst = cans[rand(cans.length)]
	cans.delete(dst)
	fg.push(dst)
	eg.push([dst, src])}
      cans += (fg + fg + Array.new(hd, src))
      src += 1
    end
    g.vertex= eg.flatten.uniq
    g.edge= eg
    g
  end

  def Ngraph.cnn(size, param=0.5)
    # connecting nearest-neighbors model introduced in Vazques 2004
    # made by nono before 050218
    # modified and added to Ngraph by fjt 2007/12/03
    vertex=[[1], [0]]
    potential=Array.new
    while(vertex.size < size) do
      if(potential.size==0) then 
	p=1.0;
      else
	p=rand();
      end
      if p < param
	### --- convert a potential edge to real
	k=rand(potential.size);
	v0,v1=potential[k];
	potential.delete_at(k);
	vertex[v0].push(v1);
	vertex[v1].push(v0);
      else
	### --- add new vertex
	jj=rand(vertex.size);
	jneighbor=vertex[jj];
	kk=vertex.size;
	vertex.push([jj]);###===vertex[kk]
	jneighbor.each do |ll|
	  potential.push([kk,ll]);
	end
	vertex[jj].push(kk);
      end
    end
    g=Ngraph.new
    g.vertex= Array.new(vertex.length){|i|i}
    i=0; g.edge= vertex.inject([]){|r, e| r+= e.map{|v|[i, v]}; i+=1; r}
    g
  end

  def Ngraph.prefo(size=100000, prob = 0.0)
    nodes=[0,1]
    links=[nodes.dup]
    (size-2).times{|i|nodes.push(i+2); vee=links.sample.sample ; links.push([vee, i+2].shuffle)
      if rand < prob
        vo=links.sample.sample
        links.push([vee, vo].shuffle) if vee != vo
      end    }
    pfn=Ngraph.new
    pfn.vertex=nodes
    pfn.edge=links
    pfn
  end

  def Ngraph.prefa(size: 100000, power: 2.0, mindeg: 1.0)
    raise "choose larger than 1 for size and 2 for power." if size < 2
    ## usage: power given by cumulative distribution function power (which is 1 smaller than the density function).
    ## any float value between 1 and 3 can be designated.


    deglb=mindeg.to_pi
    nodes = (0..deglb).to_a
    links=nodes.combination(2)

    (size-nodes.length).times{|i|
      anode=i+2
      nodes.push(anode)
      vee=nil
      mindeg.to_pi.times{vee=links.sample.sample ; links.push([vee, anode].shuffle)}

      alt = (mindeg*(2-power)/2/(power-1)).to_pi

      if alt > 0
        alt.times{
          vo=links.sample.sample
          links.push([vee, vo].shuffle) if anode != vo
        }
      else
        (-alt).times{links.delete_at(rand(links.length))}
      end
    }

    pfn=Ngraph.new
    pfn.vertex=nodes
    pfn.edge=links
    pfn
  end

  def Ngraph.complete(n)
    g = Ngraph.new
    g.vertex=Array.new(n){|i|i}
    g.edge=g.vertex.combination(2)
    g
  end

  def Ngraph.pjread(f)
    tl=File.open(f).readlines
    vn=tl.shift.split(/ +/)[1][/[0-9]+/].to_i
    nl=[]
    vn.times{nl.push(tl.shift.strip.split(/ +/)[1].gsub(/\"/, ""))}
    tl.shift
    el=tl.map{|line|line.strip.split(/ +/).map{|m|m.to_i}[0..1]}.map{|ne|[nl[ne[0]], nl[ne[1]]]}
    ret=Ngraph.new
    ret.vertex=nl
    ret.diredge=el
    ret
  end

  def Ngraph.read(f, delim = " ")
    if f.class == String
      require "pathname";
      p = Pathname.new(f).expand_path
      p.file? or raise "#{p.to_s}"
      fd = File.open(p.to_s)
    else
      raise if f.class != File
    end

    delim=Regexp.new("[" + delim + "]" + "+")
    eg = fd.readlines.map{|l|l.gsub("\n", "").split(delim)}.uniq;  fd.close
    ret = Ngraph.new
    ret.vertex = eg.flatten.uniq
    ret.diredge = eg
    ret
  end

  def initialize(obj = nil)
    @settings=eval(File.open(File.expand_path('~/.ngraph')).read)
    if obj
      i = -1
      if obj.class == Array
	@nbody = Nbody.new(obj[0].length)
	@nbody.pos = obj[i+= 1]
	@nbody.vel = obj[i+= 1]
	@nbody.acc = obj[i+= 1]
	@nbody.chg = obj[i+= 1]
	@nbody.edge = obj[i+= 1]
	@nbody.hook = obj[i+= 1]
	@nbody.mas = obj[i+= 1]
	@nbody.springlength = obj[i+= 1]
	@vertex = obj[i+= 1]
	@edge = obj[i+= 1]
	@count = obj[i+= 1]
      else
	# @nbody = Nbody.new(obj[:pos].length)
	self.vertex = obj[:vertex]
	@nbody.pos = obj[:pos]
	@nbody.vel = obj[:vel]
	@nbody.acc = obj[:acc]
	@nbody.chg = obj[:chg]
	@nbody.edge = obj[:nbedge]
	@nbody.hook = obj[:hook]
	@nbody.mas = obj[:mas]
	@nbody.airdrag = obj[:airdrag]
	@nbody.springlength = obj[:springlength]
	@edge = obj[:edge]
	@count = obj[:count]
        @mdseig = obj[:mdseig]
      end
      @nbody.tree_vector_param = @settings[:tree_vector_param]
    end
  end

  ## marshal


  def randomize(dirp=true)
    grand=Ngraph.new
    vl = self.vertex
    nl = []
    (if dirp
       [self.hailist, self.derulist].transpose.map.with_index{|e, i|[e.map{|f|f.length}, i]}
     else
       self.tonalist.map.with_index{|tl, i|[tl.length, i]}
     end).inject({}){
      |h, e|k, v = e
      if h[k]
        h[k].push(v)
      else
        h[k]=[v]
      end
      h}.to_a.transpose.last.map{|idl|
      [idl, idl.dup.shuffle].transpose.each{|p|src, dst = p; nl[dst]=vl[src]}}
    grand.vertex = nl
    grand.diredge = self.nbody.edge.map{|e|e.map{|i|nl[i]}}
    grand
  end


  def infomap(cmd:'~/bin/Infomap', opt:'-d', dir:'/tmp')
    inpath="/tmp/#{Process.pid.to_s}.net"
    outpath=inpath.gsub('.net', '.tree')
    cmd=File.expand_path(cmd)
    cmdlog=''
    ret=File.open(inpath, 'w'){|fd|
      self.infomap_net(fd)
      cmd="#{cmd} #{fd.path} #{opt} #{dir}"
      cmdlog = `#{cmd}`
      File.open(outpath).readlines[2..-1].inject([]){|ar, s|r=s.strip.split(' '); ar[r[3].to_i - 1] = [r.first.split(':')[0..-2].join(':'), r[1].to_f]; ar}
    }
    ret.unshift(cmdlog)
    File.delete(inpath)
    File.delete(outpath)
    ret
  end

  def infomap_net(out)
    head="*Vertices #{self.vertex.length}"
    boundary="*Edges #{self.edge.length}"
    out.<< head + "\n"
    out.<< self.vertex.map.with_index{|v, i|[i+1, v].join(' ')}.join("\n") +"\n"
    out.<< boundary +"\n"
    out.<< self.nbody.edge.map{|e|e.map{|i|i+1}.join(' ')}.join("\n")
    out.flush
  end

  def _dump(out)
    ## this dump order is crucial. do not touch (wara
    Marshal.dump({:pos=>@nbody.pos, :vel=>@nbody.vel,
		   :acc=>@nbody.acc,
		   :chg=>@nbody.chg,
		   :nbedge=>@nbody.edge,
		   :hook=>@nbody.hook,
		   :mas=>@nbody.mas,
		   :springlength=>@nbody.springlength,
		   :airdrag=>@nbody.airdrag,
		   :vertex=>@vertex,
		   :edge=>@edge,
		   :count=>@count,
                 :mdseig=>@mdseig}, out)
  end

  def self._load(obj)
    Ngraph.new(Marshal.load(obj))
  end

  ## data input
  def preread(path)
    require "pathname"
    File.open(Pathname.new(path).expand_path.to_s){|f|
      return f.readlines.map{|s|s.gsub("\n", "")}}
  end

  def load_vertex= (path)
    self.vertex= preread(path).uniq
    self
  end

  def load_edge= (path)
    self.edge= preread(path).map{|s|s.split(/[ ,\t]+/)}
    self
  end

  def load_diredge= (path)
    self.diredge= preread(path).map{|s|s.split(/[ ,\t]+/)}
    self
  end

  def vertex= (ary)
    @vertex= ary.uniq
    @vthash= Hash.new
    @vertex.each_with_index{|v, i|@vthash[v]= i}
    @nbody= Nbody.new(@vertex.length)
    @nbody.chg= Array.new(@vertex.length, 1/Math::sqrt(@vertex.length))
    @nbody.tree_vector_param = @settings[:tree_vector_param]
    self
  end

  def diredge=(ary) ## load edges to Nbody and setup their length, hookparams.
    @edge=ary.uniq
    # @nbody.edge=@edge.map{|e|[self.vi(e[0]), self.vi(e[1])]}
    @nbody.edge=@edge.filter{|e|flag=(i=self.vi(e[0]) and ii=self.vi(e[1])); [i, ii] if flag}
    splength=1.0/Math::log(@vertex.length)
    # splength=1.0
    # hookparam=1.0/Math::log(@vertex.length)
    hookparam=1.0
    @nbody.hook=Array.new(@edge.length, hookparam)
    @nbody.springlength=Array.new(@edge.length, splength)
    @nbody.mas=self.tonalist.map{|l|l.length.to_f + 1}
    @nbody.airdrag=Array.new(@vertex.length, 0.3)
    @tonalist=nil; tonalist
    self
  end

  alias :edge= :diredge= 

  def reinit
    @nbody = Nbody.new(@vertex.length)
    self.vertex = @vertex
    self.edge = @edge
    @count = nil
    self
  end

  def recenter
    x = @nbody.pos.map{|p|p[0]}.ave
    y = @nbody.pos.map{|p|p[1]}.ave
    z = @nbody.pos.map{|p|p[2]}.ave
    @nbody.pos = @nbody.pos.map{|p|p[0] -= x; p[1] -= y; p[2] -= z; p}
    self
  end
  
  alias :lv= :load_vertex=
  alias :le= :load_edge=

  def vertex
    @vertex
  end
  
  def edge
    @edge
  end

  def directed
    if @directed
      @directed=false
    else
      @directed=true
    end
    @directed
  end

  ## some primitive analysis

  def tonalist
    if @tonalist
      @tonalist
    else
      @tonalist=Array.new(@vertex.length){[]}
      @nbody.edge.each{|e|
	@tonalist[e[0]].push(e[1])
	@tonalist[e[1]].push(e[0])}
      @tonalist
    end
  end

  def hailist
    if @hailist
      @hailist
    else
      @hailist=Array.new(@vertex.length){[]}
      @nbody.edge.each{|e|
	@hailist[e[1]].push(e[0])}
      @hailist
    end
  end  

  def derulist
    if @derulist
      @derulist
    else
      @derulist=Array.new(@vertex.length){[]}
      @nbody.edge.each{|e|
	@derulist[e[0]].push(e[1])}
      @derulist
    end
  end

  alias :neighbours :tonalist
  
  def bfs(stt, dir=:bf, done={}, &blk)
    if stt.class != Array
      stt=[[stt]]
    else if stt.first.class != Array
           stt=[stt]
         end
    end

    stt.last.each{|e|done[e]=true}

    case dir
    when :b
      nxt=stt.last.map{|v|self.hailist[v]}.flatten.uniq.map{|e|e if not done[e]}.compact
    when :f
      nxt=stt.last.map{|v|self.derulist[v]}.flatten.uniq.map{|e|e if not done[e]}.compact
    else
      nxt=stt.last.map{|v|self.tonalist[v]}.flatten.uniq.map{|e|e if not done[e]}.compact
    end

    if nxt.length == 0 or (blk.call(stt) if blk)
      stt
    else
      bfs(stt.push(nxt), dir, done, &blk)
    end
  end

  def scc(stt)
    (self.bfs(stt, :f).flatten & self.bfs(stt, :b).flatten).uniq
  end

  def geopath(sl, dl)
    ssl=bfs(sl, -1, dl)
    ddl=bfs(dl, -1, sl).reverse
    ssl.map_with_index{|s, i|s & ddl[i]}
  end

  def geodist(src, dst, dir=:bf)
    lst=self.bfs(src, dir){|p|p.last.include?(dst)}
    lst.length if lst.last.include?(dst)
  end

  def cseg(stt, seg=nil)
    stt=[stt] if not stt.class ==Array
    seg=stt if seg.nil?
    if (nxt=(stt.map{|v|self.tonalist[v]}.flatten - seg).uniq).length == 0
      seg
    else
      cseg(nxt, seg+nxt)
    end
  end

  def connected_segments
    updatedp=true
    segvec=(0..(self.vertex.length-1)).to_a
    while updatedp
      updatedp = self.tonalist.map.with_index{|vec, i|
        min=segvec.values_at(*vec).min
        segvec[i] = min if min < segvec[i]}.inject(nil){|r,v| v or r}
    end
    segvec.map.with_index{|k, i|[k, i]}.inject({}){|h, e| k, v=e; if h[k]; h[k].push(v); else; h[k]=[v]; end; h}.to_a.transpose.last
  end

  def connected_segments_obs
    tlo=(0..self.vertex.length-1).to_a
    segments=[]
    while tlo.length > 0
      segm=cseg(tlo.first)
      segments.push(segm)
      tlo -= segm
    end
    segments
  end


  def dbfs(vl, depth=-1, vt=Array.new(@vertex.length){|i|i})
    vl=[vl] if vl.class != Array
    vt-=vl
    if vl.length == 0
      vl
    else
      if depth == 0
	[vl]
      else ## depth != 0 and some vertice given
	[vl]+dbfs(vl.map{|v|self.derulist[v]}.flatten.uniq & vt, depth-1, vt)
      end
    end
  end

  def rdbfs(vl, depth=-1, vt=Array.new(@vertex.length){|i|i})
    vl=[vl] if vl.class != Array
    vt-=vl
    if vl.length == 0
      vl
    else
      if depth == 0
	[vl]
      else ## depth != 0 and some vertice given
	[vl]+rdbfs(vl.map{|v|hailist[v]}.flatten.uniq & vt, depth-1, vt)
      end
    end
  end

  def clst(vi)
    tl = tonalist[vi] + [vi]
    tll=tl.length
    if tll> 1
      2*self.subgraph(tl).edge.length.to_f/(tll*(tll - 1))
    else
      0
    end
  end

  def smet
    tl=self.tonalist.map{|l|l.length}
    tm=tl.sort.reverse
    vn=tl.length
    en=self.edge.length
    c=smax=0; 0.upto(vn - 1){|i|
      if i > tm[i]
	lim= tm[i] - 1
      else
	lim= tm[i]
      end
      0.upto(lim){|j| smax += tm[i] * tm[j] if not i == j}}

    2 * self.nbody.edge.inject(0){|s, e| s+= (tl[e[0]] * tl[e[1]])} / smax.to_f
  end

  def average_neighbour_degree_plot
    tonalist.map{|l|
      [l.length,
	l.map{|n|tonalist[n].length}.inject(0){|r, i|r+=i}/l.length.to_f]}.\
    inject(Hash.new){|r, e|
      if r[e[0]]
	r[e[0]].push(e[1])
      else
	r[e[0]]=[e[1]]
      end
      r}.\
    to_a.\
    map{|e|
      [e[0], e[1].inject(0){|r, i|r+=i}/e[1].length.to_f]}
  end

  def adjecencymatrix
    en=self.vertex.length
    am=Array.new(en){Hash.new}
    self.nbody.edge.each{|e|am[e[0]][e[1]]=1; am[e[1]][e[0]]=1}
    am
  end

  def naconstraints(mtx= self.adjecencymatrix)  ## bogus. do not use.
    require "narray"
    en=self.vertex.length # en=g.vertex.length
    am=Array.new(en){Array.new(en, 0)}
    self.edge.each{|e|am[e[0]][e[1]] = 1} # g.edge.each{|e|am[e[0]][e[1]] = 1}
    bm=am.mplus(am.transpose); nbm = NArray.to_na(bm)
    da= (nbm * nbm).to_a
    np=bm.map{|row|row.inject(0){|v, r|r+=v}}
    sa=np 
    mx=bm.map{|row|row.max}
    m= bm.map_with_index{|row, i|row.map{|v|v/mx[i].to_f}}
    p= bm.map_with_index{|row, i|row.map{|v|v/np[i].to_f}}
    dyadredundancy=(NArray.to_na(p) * NArray.to_na(m) * NArray.to_na(da)).to_a
    r=da.mplus(dyadredundancy*(-1))
    efs=r.map{|row|row.inject(0){|v, s|s+=v}}
    efficiency=efs.map_with_index{|v, i|v/np[i]}
    pp=NArray.to_na(p) * NArray.to_na(p)
    pp2=NArray.to_na(p.mplus((pp*NArray.to_na(da)).to_a))
    (pp2*pp2).to_a.map{|row|row.inject(0){|ret, v|ret+=v}}
  end

  

  def constraints(mtx= self.adjecencymatrix)
    ata=(mtx.cmplus(mtx.cmtranspose))
    da= ata.map{|row|row.cmap{|k, v|if v > 0; 1; else; 0; end}}
    sa= ata.map{|row|ret=0; row.each{|k, v|ret+=v}; ret}
    p= ata.map_with_index{|row, i|row.cmap{|k, v|v/sa[i].to_f}}
    pp= p.cmprod(p); pp.each_with_index{|row, i|row.delete(i)}
    pp= pp.map_with_index{|row, i|row.cmap{|k, v|
        if da[i][k]
          v
        else
          0
        end}}

    #    (pp.map{|row|row.cmap{|k, v|v*v}}).map{|row|ret=0; row.each{|k, v|ret+=v}; ret}

    ((pp.cmplus(p)).map{|row|row.cmap{|k, v|v*v}})\
    .map{|row|ret=0; row.each{|k, v|ret+=v}; ret}
  end

  def tc(mtx= self.adjecencymatrix)
    ata0 =(mtx.cmplus(mtx.cmtranspose)).cm2nm
    da0= ata0.map{|row|row.map{|v|if v > 0; 1; else; 0; end}}
    np0= da0.map{|row|row.inject(0){|ret, v|ret+=v}}
    sa0= ata0.map{|row|row.inject(0){|ret, v|ret+=v}}
    p0= ata0.map_with_index{|row, i|row.map{|v|v/sa0[i].to_f}}
    pp0 = p0*p0
    ppp0= pp0.map_with_index{|row, i|row.map_with_index{|v, j|v*da0[i][j]}}
    pp20=p0.mplus(ppp0)
    (pp20.rmap{|v|v*v}).map{|row|row.inject(0){|ret, v|ret+=v}}
  end

  def lmds(plist, opt=nil, &block)
    if opt and opt[:dim]
      dim=opt[:dim]
    else
      dim=3
    end
    dmtx=plist.map{|l|
      self.bfs(l).inject_with_index(Array.new(self.vertex.length)){|v, u, i|
	u.map{|j|if block
		   v[j] = block.call(i)
		 else
		   v[j] = i
		 end}; v}}
    rowsums=dmtx.map{|row|row.inject(0){|ret, v|ret+=(v*v)}}
    colsums=dmtx.transpose.map{|row|row.inject(0){|ret, v|ret+=(v*v)}}
    tsm=rowsums.inject(0){|ret, v|ret+=v}
    en=self.vertex.length
    plen=plist.length
    bmtx=dmtx.map_with_index{|row, i|
      row.map_with_index{|v, j|
        (v*v - rowsums[i].to_f/en - colsums[j].to_f/plen + tsm.to_f/(en*plen))/(-2.0)}}
    cmtx= bmtx.map{|row|plist.map{|i|row[i]}}
    cmtx= (bmtx.mplus(bmtx.transpose)).rmap{|e|e/2}
    eig= cmtx.eigen
    pmtx= eig[0..-2].transpose
    cords=bmtx/pmtx
  end

  def dpmds(plist, opt=nil)
    if opt and opt[:dim]
      dim=opt[:dim]; p dim
    else
      dim=3
    end
    en=self.vertex.length
    eg=self.edge.length
    dmax=(Math::log(en)/(Math::log(2*eg) - Math::log(en))).ceil * 2
    mtx=plist.map{|v|
      row=Array.new(self.vertex.length, dmax)
      sp=self.bfs(v)
      goout=self.dbfs(v)
      incom=self.rdbfs(v)
      sp.each.with_index{|d, i|
        d.each{|vi|row[vi]=i}
        (goout[i]&d).each{|vi|row[vi]=2*i} if goout[i]
        (incom[i]&d).each{|vi|row[vi]=2*i} if incom[i]}
      p "#{v} done"
      row}
    p 'mtx done'
    pmds(plist, mtx, opt)
  end

  def pmds(plist, dmtx=nil, opt=nil, &block)
    if opt and opt[:dim]
      dim=opt[:dim]; p dim
    else
      dim=3
    end
    en=self.vertex.length
    eg=self.edge.length
    dmax=(Math::log(en)/(Math::log(2*eg) - Math::log(en))).ceil * 2
    if dmtx.nil?
      dmtx=plist.map{|l|
        self.bfs(l).inject_with_index(Array.new(self.vertex.length, dmax)){|v, u, i|
          u.map{|j|if block
                     v[j] = block.call(i)
                   else
                     v[j] = i
                   end}; v}}
    end
    rowsums=dmtx.map{|row|row.inject(0){|ret, v|ret+=(v*v)}}
    colsums=dmtx.transpose.map{|row|row.inject(0){|ret, v|ret+=(v*v)}}
    tsm=rowsums.inject(0){|ret, v|ret+=v}
          plen=plist.length
     bmtx=dmtx.map_with_index{|row, i|
      row.map_with_index{|v, j|
                (v*v - rowsums[i].to_f/en - colsums[j].to_f/plen + tsm.to_f/(en*plen))/(-2.0)}}
#    bt=bmtx.transpose
      b0=bmtx.prod(bmtx)
    cord=[]
    dim.times{
       e0=b0.eig_pow
       eig= Math::sqrt(Math::sqrt(e0.last))
       crd=e0.first
      proj=crd.map{|e|crd.map{|f|e*f}}.rmap{|v|v* e0.last}
      b0=b0.mplus(proj.rmap{|v|-v})
      cord.push(bmtx.transpose.prod(crd).rmap{|v|v/eig})
       }
       (3-dim).times{cord.push(Array.new(self.vertex.length, 0))} if dim < 3
    cord.transpose
  end

  def cmds(opt=nil, &block)
    if opt and opt[:dim]
      dim=opt[:dim]
    else
      dim=3
    end
    en=self.vertex.length
    eg=self.edge.length
    dmax=(Math::log(en)/(Math::log(2*eg) - Math::log(en))).ceil * 2
    dmtx=Array.new(self.vertex.length){|l|
      self.bfs(l).inject_with_index(Array.new(self.vertex.length, dmax)){|v, u, i|
	u.map{|j|if block
		   v[j] = block.call(i)
		 else
		   v[j] = i
		 end}; v}}
     sms=dmtx.map{|row|row.inject(0){|ret, v|ret+=(v*v)}}
     tsm=sms.inject(0){|ret, v|ret+=v}
    bmtx=dmtx.map_with_index{|row, i|row.map_with_index{|v, j|(v*v - sms[i].to_f/en - sms[j].to_f/en + tsm.to_f/(en*en))/(-2.0)}}
     eig=(bmtx.mplus(bmtx.transpose).rmap{|e|e/2}).eigen

      @mdseig=eig[-1]

      efecteigs=eig[-1].map_with_index{|v, i|[v, i]}.sort[-dim..-1].map{|e|e[1]}
      cord= eig[0..(self.vertex.length - 1)].transpose.map{|v|efecteigs.map{|e|v[e] * Math::sqrt(eig[-1][e])}}
      if dim < 3
        tofill = 3 - dim
        cord =  cord.map{|p|tofill.times{p.push(0.0)}; p}
      end
      if opt and opt[:setpos]
        self.pos= cord
      else
        return cord
      end
     self
  end


  def duale
    g=Ngraph.new
    g.vertex=self.nbody.edge.map{|e|e.sort}
    i=0
    r=[]
    self.tonalist.peach("edge generated"){|l|
      if l.length >= 2
	l.combination(2).each{|c|r.push([[i, c[0]].sort, [i, c[1]].sort])}
      end
      i+=1
    }
    g.edge=r
    g
  end

  def truncate ## bogus. do not use wara
    vtx=[]
    edg=[]
    0.upto(self.vertex.length - 1){|v|
      if self.tonalist[v].length > 1
	p nv = self.tonalist[v].map{|vv|[v, vv]}
	nl = nv.length
	vtx += nv
	0.upto(nl - 2){|i|  ## this is bogus. but dont know what's correct.
	  edg.push([nv[i], nv[i + 1]])}
	edg.push([nv[-1], nv[0]])
      end
    }
    self.nbody.edge.each{|ne|
      edg.push([[ne[0], ne[1]], [ne[1], ne[0]]])}
    ret=Ngraph.new; ret.vertex=vtx; ret.diredge=edg
    ret
  end

  def subgraph(verticelist)
    sg=Ngraph.new
    v=self.vertex
    if verticelist.length**2/2 > self.edge.length
      sg.vertex=verticelist.map{|i|v[i]}
      vh=sg.vertex.inject({}){|h,v|h[v]=true;h}
      sg.edge=self.edge.filter{|e|e if vh[e.first] and vh[e.last]}
    else
      ll=self.derulist
      nv=[]
      ne=[]
      verticelist.each{|i|nv.push(v[i]); ne+=((ll[i] & verticelist).map{|j|[v[i], v[j]]})}
      sg.vertex=nv; sg.diredge=ne
    end
    sg
  end

  def update(update)
    gu=Ngraph.new
    vertices=update[:vertice].uniq
    edge=update[:edge].uniq
    gu.vertex=vertices
    edge=edge.filter{|e|e if gu.vi(e.first) and gu.vi(e.last)}
    gu.diredge=edge
    found=Array.new(vertices.length){nil}
    sp=self.pos
    gup=gu.pos
    vertices.each.with_index{|v,i|
      if j=self.vi(v)
        gup[i]=sp[j]; found[i]=true
      end}

    gut=gu.tonalist
    cseg=found.map.with_index{|e, i|i if e}.compact
    gu.bfs(cseg).each_cons(2){|slice|
      slice.last.each{|vi|
        cp=gut[vi].filter{|i|i if found[i]}.map{|i|gup[i]}.transpose.map{|v|v.ave}
        gup[vi]=cp
        found[vi]=true
      }}
    gu.pos=gup
    gu
  end


  # def subgraph(vertice, opt={})
  #   g = Ngraph.new
  #   neg = [];
  #   tl = self.tonalist
  #   ne=self.nbody.edge
  #   vv=vertice.dup
  #   # vertice.each{|i|neg+= ((tl[i] & vertice).map{|v|[i, v]})}
  #   vertice.each{|i|ne = (tl[i] & vv)
  #     if opt[:mst]
  #       vv-=tl[i]
  #     else
  #       vv-= [i]
  #     end
  #     neg+= ne.map{|v|[i, v]}}


  #   neg = neg.uniq
  #   nv = vertice
    
  #   if opt[:name]
  #     v = self.vertex; nv = vertice.map{|i|v[i]}
  #     neg = neg.rmap{|i|v[i]}
  #   end

  #   g.vertex = nv

  #   g.edge = neg
  #   p = self.pos; g.pos = vertice.map{|i|p[i]}
  #   g
  # end

  def connected?
    self.bfs(0).flatten.length == self.vertex.length
  end

  def add(vertices, edges)
    
  end

  def subt(vindices, edges)
  end

  def +(graph)
    if graph.class != Ngraph
      raise
    else
      g=Ngraph.new
      g.vertex=(self.vertex + graph.vertex).uniq
      g.edge=(self.edge + graph.edge).uniq if self.edge and graph.edge
      g
    end
  end

  def inter(n)
    g=Ngraph.new
    lv=self.vertex.length
    v=Array.new(self.vertex.length){|i|i}
    eg=[]
    i=0; 
    self.nbody.edge.peach{|e|
      ary=[e[0]]+Array.new(n){|j|lv+n*i+j}+[e[1]]
      v+=Array.new(n){|j|lv+n*i+j}
      (n+1).times{|k|eg.push([ary[k], ary[k+1]])}
      i+=1
    }
    g.vertex=v; g.edge=eg
    g
  end

  
  ## interface functions
  def frog(dt, &block)
    if block
      block.call(@nbody)
    end
    @nbody.frog(dt.to_f)
    count(dt)
    self
  end

  def pure_frog(dt, &block)
    if block
      block.call(@nbody)
    end
    @nbody.pure_frog(dt.to_f)
    count(dt)
    self
  end

  def mgn(value=1.0, dir=[0.0, -1.0, 0.0])
    @nbody.mgn(value.to_f, dir.rmap{|f|f.to_f})
    self
  end

  def mgn_ruby(value, dir)
    ## dir=[0, -1, 0] makes incomming node go down; value=0.1 not bad.
    dir.class == Array or raise "give direction with array"
    dim = @nbody.pos[0].length
    dir = dir.map{|d|d/dir.abs}
    ac = @nbody.acc; pos = @nbody.pos; mas = @nbody.mas
    @nbody.edge.each{|e|
      m0 = mas[e[0]]
      m1 = mas[e[1]]
      ev = Array.new(dim){|i| pos[e[1]][i] - pos[e[0]][i]}
      if ev.abs > 0
	ev = ev.map{|eg|eg/ev.abs}
	rv = Array.new(dim){|i| dir[i] - ev[i]}
	ac[e[0]] = Array.new(dim){|i| ac[e[0]][i] - (value * rv[i])/m0}
	ac[e[1]] = Array.new(dim){|i| ac[e[1]][i] + (value * rv[i])/m1}
      end
      }
    @nbody.acc = ac
    self
  end

  def distance(v1, v2)
    dist(vi(v1), vi(v2))
  end
  
  def dist(vi1, vi2)
    dim = @nbody.pos[0].length; d = 0; pos = @nbody.pos
    0.upto(dim - 1){|i|
      d +=  (pos[vi1][i] - pos[vi2][i])**2
    }
    Math::sqrt(d)
  end

  def epotential
    mass = self.mas
    self.acc.map_with_index{|a, i|Math::sqrt(a.abs * mass[i])}.inject(0){|r, p|r += p}
  end

  def bpotential
    sl = self.springlength
    hk = self.hook
    pss = self.pos
    @nbody.edge.inject_with_index(0){|ret, eg, i|
      ret += hk[i] * \
      ((pss[eg[0]].map_with_index{|e, i|e - pss[eg[1]][i]}.abs - sl[i]) ** 2) / 2
    }
  end

  def kinetique
    mass = self.mas
    self.vel.map_with_index{|v, i|mass[i] * (v.abs ** 2) / 2.0}.inject(0){|k, i|k+= i}
  end

  def pot
    self.epotential + self.bpotential
  end

#   def fr(dt, &blk)
#     ps = @nbody.pos
#     ac = @nbody.acc
#     0.upto(@vertex.length - 1){|i|
#       a = ac[i].abs
#       0.upto(2){|j|
# 	ps[i][j] += blk.call(ac[i][j], a)*dt}}
#     @nbody.pos = ps
#     count(dt)
#     self
#   end

  def fr(dt)
    ps = @nbody.pos
    ac = @nbody.acc
    @nbody.pos = ps.map_with_index{|p, i|
      a=ac[i]
      b=a.abs
      p.map.with_index{|pp, ii|
	pp += dt * a[ii] / b}}
    count(dt)
    self
  end

  def converge_pos(dt = 0.1, ratio = 1.1, &block)
    flag = true
    count = Array.new
    pp = self.pos; pt = self.pot

    while flag
      block.call(dt)
      pc = self.pos
      flag = pc.map_with_index{|e, i|
	e.map_with_index{|ee, ii|
	  pp[i][ii] - ee}.abs}.max > ( pc.map{|p|p.abs}.ave / 10000)
      count.push(self.pot < pt); count.shift if count.length > 3
      if count.inject(true){|r, e|r and e}
	dt = dt * ratio
      else
	dt = dt / ratio
      end
      
      pp = self.pos; pt = self.pot
    end
    self
  end

  def pots
    if @pots
      @pots
    else
      @pots=Array.new(5, 0)
    end
  end      

  def count(dt = nil)
    if dt
      if @count
	@count[0]+= 1
	@count[1]+= dt
      else
	@count = [1, dt]
      end
    end
    @count
  end


  ["pos",
    "vel",
    "acc",
    "mas",
    "chg",
    "hook",
    "springlength",
    "airdrag"].each{|mname|
    define_method(mname){@nbody.send(mname)}
    mname2= mname+"="
    define_method(mname2){|arg|
      case mname
      when "hook", "springlength"
	len = @edge.length
      else
	len = @vertex.length
      end
      if arg.class != Array
	arg = Array.new(len, arg)
      end      
      @nbody.send(mname2, arg.rmap{|i|i.to_f})}}

  alias :spln :springlength
  alias :spln=  :springlength=
  alias :air :airdrag
  alias :air= :airdrag=
  
  
  def force
    @nbody.coulomb_force
    self
  end

  def bane
    @nbody.bane
    self
  end

  def nbody
    @nbody
  end

  def vi(v)
    @vthash[v]
  end

  def pos2vrml
    diam = 0.01
    "#VRML V1.0 ascii\nSeparator { Switch {\n" + \
    self.pos.map{|ps|
      [ps,
	[ps[0] + diam, ps[1], ps[2]],
	[ps[0], ps[1] + diam, ps[2]],
	[ps[0], ps[1], ps[2] + diam],
	[ps[0] + diam, ps[1] + diam, ps[2]],
	[ps[0] + diam, ps[1], ps[2] + diam],
	[ps[0], ps[1] + diam, ps[2] + diam],
	ps.map{|pos|pos + diam}]
    }.map{|cube|cube.map{|v|v.join(" ")}.join(",\n") + ",\n"}.map_with_index{|cc, i|
      "DEF Vertex#{i} Separator { Coodinate3 { point [ #{cc} ] } IndexedFaceSet {coordIndex [ 0, 1, 2, 3, -1,
				 4, 7, 6, 5, -1,
				 0, 4, 5, 1, -1,
				 1, 5, 6, 2, -1,
				 2, 6, 7, 3, -1,
				 4, 0, 3, 7, -1,]}}\n"}.join("\n") + \
    "}"
  end


  def plain_view(host="localhost", port=9876)
    push(host, port){
      out=""
      out+=@nbody.pos.length.to_s + "\n"
      out+=@nbody.pos.map{|pos|pos.join(" ")}.join("\n")+"\n"
      out+=@nbody.edge.map{|e|e.join(" ")}.join("\n")+"\n"
      out}
  end

  def colour(&colour)
    @colour=colour
  end

  def size(&size)
    @size=size
  end

  def label(&label)
    @label=label
  end

  def drawc(arg={})
    require "cairo"
    size=(arg[:size] or 2000)
    path=(arg[:path] or "./hoge.png")
    edge=(arg[:edge] or nil)
    scale=(arg[:scale] or 1.0)
    sp=self.pos
    if arg[:vertice].nil?
      poss= sp
      vertice =[sp.map{|p|p.last}, Array.new(self.vertex.length){|i|i}].transpose.sort.transpose.last
    else
      vertice=arg[:vertice].compact
      vertice =vertice.map{|i|[sp[i].last, i]}.sort.transpose.last
      vh=vertice.inject({}){|h, e|h[e]=true; h}
    end ## sort vertice indices along z axis position

    sgm=Math::sqrt(sp.map{|p|p.inject(0){|a,b|a+b*b}}.ave)
    conv=lambda{|p|p*size/(5.0*sgm/scale)+size/2.0}
    poss= sp.rmap{|v|conv.call(v)}

    Cairo::ImageSurface.new(size, size){|surface|
      cntxt=Cairo::Context.new(surface)
      cntxt.set_source_rgb(0, 0, 0)
      cntxt.paint

      if edge
        cntxt.set_source_rgba(1.0, 1.0, 1.0, 0.2) ## opaque white is (1, 1, 1, 1)

        if arg[:vertice].nil?
          elist= self.nbody.edge
        else
          elist= self.nbody.edge.filter{|e|e if vh[e.first] and vh[e.last]}
        end
        elist.each{|e|
          cntxt.move_to(*(poss[e[0]].first(2)))
          cntxt.line_to(*(poss[e[1]].first(2)))
          cntxt.set_line_width(1.0)
          cntxt.stroke
        }
      end

      cntxt.set_source_rgba(1, 0, 0, 1.0)
      vertice.each{|vi|
        gp= poss[vi]
        ecs, wai=gp[0..1]
        if @size
          diam=@size.call(vi)
        else
          diam= (size*size)/ (self.vertex.length * 500.0)
        end
        diam= 1.0 if diam < 1.0
        if @colour
          rgb= @colour.call(vi)
          if rgb.length == 3
            rgb.push(1.0) ## always with alpha channel
          end
          cntxt.set_source_rgba(*rgb)
        end
        cntxt.circle(ecs, wai, diam).fill
        if @label
          str= @label.call(vi).to_s
          fsize = if diam < 10 ; 10 ;else; diam ;end
          cntxt.set_font_size(fsize)
          cntxt.move_to(ecs+diam, wai)
          cntxt.show_text(str)
        end
      }

      surface.write_to_png(path)
    }
  end

  def draw2d (arg={})
    size=(arg[:size] or 2000)
    path=(arg[:path] or "/tmp/hoge")
    edge=(arg[:edge] or nil)
    scale=(arg[:scale] or 1.0)
    require "GD"
    sp=self.pos
    poss=sp.dup
    if arg[:vertice].nil?
#      vertice=Array.new(self.vertex.length){|i|i}
      poss= sp
      vertice =[sp.map{|p|p.last}, Array.new(self.vertex.length){|i|i}].transpose.sort.transpose.last
    else
      vertice=arg[:vertice].compact
      vertice =vertice.map{|i|[sp[i].last, i]}.sort.transpose.last
      vh=vertice.inject({}){|h, e|h[e]=true; h}
    end
    # if opts[:scale] == :local
    #   gw= poss.transpose[0].max - poss.transpose[0].min
    #   gh= poss.transpose[1].max - poss.transpose[1].min
    #   xs= poss.transpose[0].min; ys= poss.transpose[1].min
    # else
    #   gw= sp.transpose[0].max - sp.transpose[0].min
    #   gh= sp.transpose[1].max - sp.transpose[1].min
    #   xs= sp.transpose[0].min; ys= sp.transpose[1].min
    # end
    # poss= poss.map{|p|p[0]-=xs; p[1]-=ys; p}

    # if (gw/gh) > (width.to_f / height) ## ratio to width, else height
    #   conv= lambda{|p|0.9*(p*width / gw) + width/20.0}
    # else
    #   conv= lambda{|p|0.9*(p*height / gh) + height/20.0}
    # end
    p scale
    p sgm=Math::sqrt(sp.map{|p|p.inject(0){|a,b|a+b*b}}.ave)
    conv=lambda{|p|p*size/(5.0*sgm/scale)+size/2.0}
    poss= sp.rmap{|v|conv.call(v)}

    @drawing= GD::Image.newTrueColor(size, size)
    @drawing.colorAllocate(0, 0, 0)
    red= @drawing.colorAllocate(255, 0, 0)
    gray= @drawing.colorAllocateAlpha(255,255,255,112)
    if edge
      if arg[:vertice].nil?
        elist= self.nbody.edge
      else
        elist= self.nbody.edge.filter{|e|e if vh[e.first] and vh[e.last]}
      end
      elist.each{|e|
        @drawing.line(poss[e[0]][0], poss[e[0]][1], poss[e[1]][0], poss[e[1]][1], gray)}
    end

#    vertice.each_with_index{|vi, i|
    vertice.each{|vi|
      gp= poss[vi]
      if @size
        diam=@size.call(vi)
      else
        # diam= size/(5*Math::sqrt(self.vertex.length))
        diam= (size*size)/ (self.vertex.length * 500.0)
      end
      diam= 1.0 if diam < 1.0
      if @colour
        rgb= @colour.call(vi).map{|v|v.to_i}
 #        colour= @drawing.colorAllocate(rgb[0], rgb[1], rgb[2])
        if rgb.length > 3
          colour= @drawing.colorAllocateAlpha(*rgb)
        else
          colour= @drawing.colorAllocate(*rgb)
        end
      else
        colour= red
      end
      @drawing.filledEllipse(gp[0], gp[1], diam, diam, 0, 360, colour)
      if @label
        str= @label.call(vi)
        strcol=@drawing.colorAllocate(255, 255, 255)
        @drawing.stringTTF(strcol, "/usr/share/fonts/truetype/kochi/kochi-gothic.ttf", diam, 0, gp[0]+diam, gp[1], str)
      end
}

    File.open(path, "w"){|out|@drawing.png(out)}
    @drawing.destroy
  end

  def dumpview
    u=`whoami`[0..-2]
    path="/tmp/" + u + Time.now.strftime("%s").to_s
    File.open(path, "w"){|o|Marshal.dump(self, o)}
    system("~/bin/grv #{path} #{$KCODE}")
    File.open(path){|fd|
      g=Marshal.load(fd.read)
      self.pos = g.pos
      self.vel = g.vel
      self.acc = g.acc
      @count = g.count}
    require "fileutils"; FileUtils.rm(path)
  end
  
  def view_with_size_colour(host="localhost", port=9876)
    if @size.nil?
      size=Array.new(@vertex.length, 1.0)
    else
      size=Array.new(@vertex.length){|i|@size.call(i)}
    end
    if @colour.nil?
      colour=Array.new(@vertex.length, [1.0, 1.0, 0.0])
    else
      colour=Array.new(@vertex.length){|i|@colour.call(i)}
    end
    push(host, port){
      print "\n"
      out=""
      out+=@nbody.pos.length.to_s + "\n"
      @nbody.pos.each_with_index{|pos, i|out+=pos.join(" ")+" "+
	  size[i].to_s+" "+
	  colour[i].join(" ")+"\n"
	print i+1, "/", @vertex.length, "\r"
      }
      out+="\n"
      out+=@nbody.edge.map{|e|e.join(" ")}.join("\n")+"\n"
      out}
  end

  alias :view :plain_view
  alias :nview :view_with_size_colour

  def push(host, port, &block)
    require "socket"
    user=`whoami`
    soc=TCPSocket.open(host, port)
    soc.write("user #{user}\n")
    soc.write("put\n")
    soc.gets
    case $_
    when /^ok/i
      out=block.call
      outsize=out.length.to_s
      soc.write(outsize+"\n")
      soc.write(out)
      soc.gets
      if $_[/ok/i]
	soc.close
      else
	raise
      end
    else
      print "NG\n"
      soc.close
    end
    true
  end
end
