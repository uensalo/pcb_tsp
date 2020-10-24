// Author: Unsal Ozturk
int MODE = 1; // 0 for Manhattan Distance, 1 for Floyd-Warshall

class Point {
  int x;
  int y;
  Point(int x, int y) {
    this.x = x;
    this.y = y;
  }
  
  @Override
  public String toString() {
    return "Point("+x+","+y+")";
  }
  
  @Override
  public boolean equals(Object other){
    if(other instanceof Point) {
      return ((Point)other).x == x && ((Point)other).y == y;
    }
    return false;
  }
}

class Block {
  Point aa;
  Point bb;
  Block(Point aa, Point bb) {
    this.aa = aa;
    this.bb = bb;
  }
}

class Path {
  ArrayList<Point> vertices;
  Path() {
    this.vertices = new ArrayList<Point>();
  }
  Path(ArrayList<Point> vertices) {
    this.vertices = new ArrayList<Point>();
    for(Point pt: vertices) {
      Point pnew = new Point(pt.x, pt.y);
      this.vertices.add(pnew);
    }
  }
  @Override
  public String toString() {
    String str = "";
    for(Point v: vertices) {
      str += v.toString() + "->";
    }
    return str;
  }
  public int length() {
    return vertices.size();
  }
}

class Vertex {
  ArrayList<Integer> adjlist;
  Vertex() {
    adjlist = new ArrayList<Integer>();
  }
}

class Graph {
  Vertex[] vertices;
  Graph(int N) {
    vertices = new Vertex[N];
    for(int i = 0; i < N; i++) {
      vertices[i] = new Vertex();
    }
  }
  void addEdge(int u, int v) {
    vertices[u].adjlist.add(v);
    vertices[v].adjlist.add(u);
  }
  
  Path getCycle(int src) {  
    ArrayList<Integer> s = new ArrayList<Integer>();
    int dst = vertices[src].adjlist.get(0);
    vertices[src].adjlist.remove(0);
    vertices[dst].adjlist.remove(Integer.valueOf(src));
    s.add(src);
    boolean[] visited = new boolean[vertices.length];
    int[] parent = new int[vertices.length];
    parent[src] = -1;
    while(!s.isEmpty()) {
      int cur = s.get(s.size() - 1);
      s.remove(s.size() - 1);
      if(!visited[cur]) {
        visited[cur] = true;
      }
      for(int u : vertices[cur].adjlist) {
        if(!visited[u]) {
          s.add(u);
          parent[u] = cur;
        }
      }
    }
    
    Path cycle = new Path();
    int cur = dst;
    while(parent[cur] != -1) { 
      Point p = new Point(parent[cur] + 1,cur + 1);
      cycle.vertices.add(0, p);
      cur = parent[cur];
    }
    vertices[src].adjlist.add(0, dst);
    vertices[dst].adjlist.add(src);
    cycle.vertices.add(new Point(dst + 1, src + 1));
    return cycle;
  }
} 

final int w = 800;
final int h = 800;
String PATH_TABLE_FW = "table_fw_paths.csv";
String PATH_TABLE_TC = "table_tc_paths.csv";
String PATH_SOLN_FW = "tsp_fw_soln.csv";
String PATH_SOLN_TC = "tsp_tc_soln.csv";
String txtedges;
int nPoints;
int nEdges;
int nRows;
int nCols;
int[] solution;
int curSolnVal;
Table tablePoints;
Table tableBlocks;
Table tablePaths;
Table tableSoln;
ArrayList<Point> points;
ArrayList<Block> blocks;
ArrayList<Path> paths;
ArrayList<Integer> solutionEdgeIdx;
ArrayList<Integer> sortedSolutionEdgeIdx;
ArrayList<Integer> sortedNodeIndices;
int curEdge = 0;

void init() {
  tablePoints = loadTable("points.csv");
  curSolnVal = 0;
  int t = 0;
  nRows = -1;
  nCols = -1;
  points = new ArrayList<Point>();
  blocks = new ArrayList<Block>();
  paths = new ArrayList<Path>();
  sortedNodeIndices = new ArrayList<Integer>();
  for(TableRow row : tablePoints.rows()) {
    if(t == 0) { t++; continue;}
    int x = row.getInt(0);
    int y = row.getInt(1);
    if (x > nRows) nRows = x;
    if (y > nCols) nCols = y;
    Point p = new Point(x - 1,y - 1);
    points.add(p);
    t++;
  }
  nPoints = points.size();
  solution = new int[points.size() * (points.size() - 1) / 2];
  
  tableBlocks = loadTable("rects.csv");
  t = 0;
  for(TableRow row : tableBlocks.rows()) {
    if(t == 0) { t++; continue;}
    Point aa = new Point(row.getInt(0) - 1, row.getInt(1) - 1);
    Point bb = new Point(row.getInt(2) - 1, row.getInt(3) - 1);
    blocks.add(new Block(aa,bb));
    t++;
  }
  
  tablePaths = MODE == 0 ? loadTable(PATH_TABLE_TC) : loadTable(PATH_TABLE_FW);
  t = 0;
  for(TableRow row : tablePaths.rows()) {
    if(t == 0) { t++; continue; }
    int noVertices = 0;
    ArrayList<Point> vertices = new ArrayList<Point>();
    for(int i = 0; i < row.getColumnCount(); i++) {
      if(row.getInt(i) == 0) {
        break;
      }
      noVertices++;
    }
    noVertices /= 2;
    for(int i = 0; i < noVertices; i++) {
      Point p = new Point(row.getInt(i) - 1, row.getInt(i+noVertices) - 1);
      vertices.add(p);
    }
    Path path = new Path(vertices);
    paths.add(path);
    t++;
  }
  
  tableSoln = MODE == 0 ? loadTable(PATH_SOLN_TC) : loadTable(PATH_SOLN_FW);
  t = 0;
  solutionEdgeIdx = new ArrayList<Integer>();
  for(TableRow row : tableSoln.rows()) {
    solution[t] = row.getInt(0);
    if(row.getInt(0) == 1) {
      solutionEdgeIdx.add(t);
    }
    t++;
  }
  
  //generate n choose k pairs
  ArrayList<Point> pair_list = new ArrayList<Point>();
  for(int i = 1; i <= nPoints - 1; i++) {
    for(int j = i + 1; j <= nPoints; j++) {
      pair_list.add(new Point(i,j));
    }
  }
  
  Graph g = new Graph(nPoints);
  ArrayList<Point> pairs_in_solution = new ArrayList<Point>();
  for(int i = 0; i < solutionEdgeIdx.size(); i++) {
    pairs_in_solution.add(pair_list.get(solutionEdgeIdx.get(i)));
    g.addEdge(pairs_in_solution.get(i).x - 1,pairs_in_solution.get(i).y - 1);
  }
  
  Path pairs_sorted_by_cycle = g.getCycle(0);
  sortedSolutionEdgeIdx = new ArrayList<Integer>();
  txtedges = "";
  for(Point p : pairs_sorted_by_cycle.vertices) {
    int max = p.x > p.y ? p.x : p.y;
    int min = p.x > p.y ? p.y : p.x;
    print(p.x + " - ");
    sortedNodeIndices.add(p.x);
    Point augP = new Point(min,max);
    int idx = pair_list.indexOf(augP);
    sortedSolutionEdgeIdx.add(idx);
  }
  sortedNodeIndices.add(1);
  println(1);
}

void setup() {
  size(800,800);
  init();
}  

void draw() {
  background(255);
  drawGrid();
  fill(0);
  String mode = "Paths in between holes are allowed only through the shortest Manhattan Distance";
  if(MODE != 0)
    mode = "Paths in between holes are always allowed, regardless of the shortest Manhattan Distance";
  text("Current Mode: " + mode, 50,20);
  text("Current Path Length: " + curSolnVal, 50,40);
  text("Current Number of Edges: " + curEdge, 50,60);
  text("Edges: " + txtedges, 50,80); 
  text("LEFT/RIGHT arrow key to DECREASE/INCREASE number of edges by one.", 50,730);
  text("UP arrow key to change drill movement mode.", 50,750);
}

void keyPressed() {
  if (keyCode == RIGHT) {
    curEdge = ++curEdge > 50 ? 0 : curEdge;
  } else if (keyCode == LEFT) {
    curEdge = --curEdge < 0 ? 50 : curEdge;
  } else if (keyCode == UP || keyCode == DOWN) {   
    MODE = MODE == 0 ? 1 : 0;
    println("Changing mode to " + MODE + "...");
    init();
  } 
}

void drawGrid() {
  float startX = 100;
  float startY = 120;
  float endX = 700;
  float endY = 720;
  
  float dr = (endY - startY) / nRows;
  float dc = (endX - startX) / nCols;
  
  fill(255);
  for(int i = 0; i < nRows; i++) {
    for(int j = 0; j < nCols; j++) {
      circle(startX + j * dc, startY + i * dr, 10);
    }
  }
  
  fill(204,102, 0);
  for(Block b: blocks) {
    Point aa = b.aa;
    Point bb = b.bb;
    for (int i = aa.x; i <= bb.x; i++) {
      for (int j = aa.y; j <= bb.y; j++) {
        circle(startX + j * dc, startY + i * dr, 10);
      }
    }
  }
  
  txtedges = "";
  fill(102, 102, 0);
  int solval = 0;
  for(int i = 0; i < curEdge; i++) {
    if(i != 0 && i % 28 == 0) {
      txtedges = txtedges + "\n";
    }
    txtedges = txtedges + sortedNodeIndices.get(i) + " - ";
    Path path = paths.get(sortedSolutionEdgeIdx.get(i));
      for(Point p : path.vertices) {
        solval++;
        circle(startX + p.y * dc, startY + p.x * dr, 10);
      }
      solval--;
  }
  txtedges = txtedges + (sortedNodeIndices.get(curEdge));
  curSolnVal = solval;

  int t = 0;
  for(Point p: points) {
    if(t == 0) {t++; fill(0,0, 102);}
    else fill(0,102, 0);
    circle(startX + p.y * dc, startY + p.x * dr, 10);
    t++;
  }
}
