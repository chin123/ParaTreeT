#ifndef SIMPLE_TREEPIECE_H_
#define SIMPLE_TREEPIECE_H_

#include "simple.decl.h"
#include "common.h"
#include "templates.h"
#include "ParticleMsg.h"
#include "MultiMsg.h"
#include "Node.h"
#include "Utility.h"
#include "Reader.h"
#include "CentroidVisitor.h"
#include "GravityVisitor.h"

#include <queue>
#include <set>
#include <stack>
#include <queue>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

extern CProxy_Reader readers;
extern int max_particles_per_leaf;
extern int decomp_type;
extern int tree_type;
extern CProxy_TreeElement<CentroidData> centroid_calculator;
extern CProxy_CacheManager<CentroidData> centroid_cache;
extern CProxy_Main mainProxy;

template <typename Data>
class TreePiece : public CBase_TreePiece<Data> {
  std::vector<Particle> particles;
  std::vector<Node<Data>*> leaves;
  int n_total_particles;
  int n_treepieces;
  int particle_index;
  int n_expected;
  Key tp_key; // should be a prefix of all particle keys underneath this node
  Node<Data>* root;
  Node<Data>* root_from_tp_key;
  std::vector<std::set<Key>> curr_nodes;
  std::vector<Node<Data>*> trav_tops;
  CProxy_TreeElement<Data> global_data;
  CProxy_CacheManager<Data> cache_manager;
  std::set<Key> curr_waiting;
  int num_done;
  // debug
  std::vector<Particle> flushed_particles;

public:
  TreePiece(const CkCallback&, int, int); 
  void receive(ParticleMsg*);
  void check(const CkCallback&);
  void triggerRequest();
  void build(const CkCallback&);
  bool recursiveBuild(Node<Data>*, bool);
  template<typename Visitor>
  void upOnly(CProxy_TreeElement<Data>);
  template<typename Visitor>
  void startDown(CProxy_CacheManager<Data>);
  template<typename Visitor>
  void startUpAndDown(CProxy_CacheManager<Data>);
  template<typename Visitor>
  void requestNodes(Key, CProxy_CacheManager<Data>, int);
  template<typename Visitor>
  void goDown(Key); 
  void print(Node<Data>*);
  Node<Data>* findNode(Key, bool start_from_tp_key = false);
  void perturb (Real timestep);
  void flush(CProxy_Reader);
  void rebuild(const CkCallback&);

  // debug
  void checkParticlesChanged(const CkCallback& cb) {
    bool result = true;
    if (particles.size() != flushed_particles.size()) {
      result = false;
      this->contribute(sizeof(bool), &result, CkReduction::logical_and_bool, cb);
      return;
    }
    for (int i = 0; i < particles.size(); i++) {
      if (!(particles[i] == flushed_particles[i])) {
        result = false;
        break;
      }
    }
    this->contribute(sizeof(bool), &result, CkReduction::logical_and_bool, cb);
  }
};

template <typename Data>
TreePiece<Data>::TreePiece(const CkCallback& cb, int n_total_particles_, int n_treepieces_) : n_total_particles(n_total_particles_), n_treepieces(n_treepieces_), particle_index(0) {
  if (decomp_type == OCT_DECOMP) {
    // OCT decomposition
    n_expected = readers.ckLocalBranch()->splitters[this->thisIndex].n_particles;
    tp_key = readers.ckLocalBranch()->splitters[this->thisIndex].tp_key;
  }
  else if (decomp_type == SFC_DECOMP) {
    // SFC decomposition
    n_expected = n_total_particles / n_treepieces;
    if (this->thisIndex < (n_total_particles % n_treepieces))
      n_expected++;
      // TODO tp_key needs to be found in local tree build
  }
  this->contribute(cb);
}

template <typename Data>
void TreePiece<Data>::receive(ParticleMsg* msg) {
  // copy particles to local vector
  particles.resize(particle_index + msg->n_particles);
  std::memcpy(&particles[particle_index], msg->particles, msg->n_particles * sizeof(Particle));
  particle_index += msg->n_particles;
  delete msg;
}

template <typename Data>
void TreePiece<Data>::check(const CkCallback& cb) {
  if (n_expected != particles.size()) {
    CkPrintf("[TP %d] ERROR! Only %d particles out of %d received\n", this->thisIndex, particles.size(), n_expected);
    CkAbort("Failure on receiving particles");
  }
  this->contribute(cb);
}

template <typename Data>
void TreePiece<Data>::triggerRequest() {
  readers.ckLocalBranch()->request(this->thisProxy, this->thisIndex, n_expected);
}

template <typename Data>
void TreePiece<Data>::build(const CkCallback &cb){
  // sort particles received from readers
  std::sort(particles.begin(), particles.end());

  // create global root and recurse
#ifdef DEBUG
  CkPrintf("[TP %d] key: 0x%" PRIx64 " particles: %d\n", this->thisIndex, tp_key, particles.size());
#endif
  if (decomp_type == OCT_DECOMP) {
    root = new Node<Data>(1, 0, particles.size(), &particles[0], 0, n_treepieces - 1, NULL);
    recursiveBuild(root, false);
  } else if (decomp_type == SFC_DECOMP) {
    root = new Node<Data>(1, 0, particles.size(), &particles[0], 0, n_treepieces - 1, NULL);
    root->type = Node<Data>::Shared;

    std::queue<Node<Data>*> queue;
    queue.push(root);

    while (!queue.empty()) {
      Node<Data>* node = queue.front();
      queue.pop();

      // Create each child
      for (Key child_key = (node->key >> LOG_BRANCH_FACTOR);
           child_key < (node->key >> LOG_BRANCH_FACTOR) + BRANCH_FACTOR;
           child_key++)
      {
        // Decide node type
        int start_key = Utility::removeLeadingZeros(child_key);
        int end_key = Utility::removeLeadingZeros(child_key) + 1;
        int start_idx = Utility::binarySearchG(start_key, node->particles, 0, node->n_particles);
        int end_idx = Utility::binarySearchGE(end_key, node->particles, 0, node->n_particles);
        
        Node<Data>* child = new Node<Data>(
          child_key, node->depth+1, end_idx-start_idx+1,
          node->particles+start_idx-1, 0, 0, node
        );
        if (end_idx <= 0 || start_idx >= node->n_particles) {
          child->type = Node<Data>::Remote;
        } else if (start_idx > 0 && end_idx < node->n_particles) {
          child->type = Node<Data>::Local;
        } else {
          child->type = Node<Data>::Shared;
        }
        
        node->children.push_back(child);

        bool is_light = (node->n_particles <= ceil(BUCKET_TOLERANCE * max_particles_per_leaf));
        if (node->type == Node<Data>::Shared || (node->type == Node<Data>::Local && is_light)) {
          queue.push(child);
        }
      }

    }
  }

  this->contribute(cb);
}

template <typename Data>
bool TreePiece<Data>::recursiveBuild(Node<Data>* node, bool saw_tp_key) {
#if DEBUG
  //CkPrintf("[Level %d] created node 0x%" PRIx64 " with %d particles\n",
    //  node->depth, node->key, node->n_particles);
#endif
  // store reference to splitters
  //static std::vector<Splitter>& splitters = readers.ckLocalBranch()->splitters;

  if (tree_type == OCT_TREE) {
    // check if we are inside the subtree rooted at the treepiece's key
    if (!saw_tp_key) {
      saw_tp_key = (node->key == tp_key);
    }

    bool is_light = (node->n_particles <= ceil(BUCKET_TOLERANCE * max_particles_per_leaf));
    bool is_prefix = Utility::isPrefix(node->key, tp_key);
    /*
    int owner_start = node->owner_tp_start;
    int owner_end = node->owner_tp_end;
    bool single_owner = (owner_start == owner_end);
    */

    if (saw_tp_key) {
      // XXX if we are under the TP key, we can stop going deeper if node is light
      if (is_light) {
        if (node->n_particles == 0)
          node->type = Node<Data>::EmptyLeaf;
        else {
          node->type = Node<Data>::Leaf;
          leaves.push_back(node);
        }
        return true;
      }
    }
    else { // still on the way to TP, we could have diverged
      if (!is_prefix) {
        // diverged, should be remote
        node->type = Node<Data>::Remote;

        return false;
      }
    }

    // for all other cases, we go deeper

    /*
    if (is_light) {
      if (saw_tp_key) {
        // we can make the node a local leaf
        if (node->n_particles == 0)
          node->type = Node<Data>::EmptyLeaf;
        else
          node->type = Node<Data>::Leaf;

          return true;
      }
      else if (!is_prefix) {
        // we have diverged from the path to the subtree rooted at the treepiece's key
        // so designate as remote
        node->type = Node<Data>::Remote;

        CkAssert(node->n_particles == 0);
        CkAssert(node->n_children == 0);

        if (single_owner) {
          int assigned = splitters[owner_start].n_particles;
          if (assigned == 0) {
            node->type = Node<Data>::RemoteEmptyLeaf;
          }
          else if (assigned <= BUCKET_TOLERANCE * max_particles_per_leaf) {
            node->type = Node<Data>::RemoteLeaf;
          }
          else {
            node->type = Node<Data>::Remote;
            node->n_children = BRANCH_FACTOR;
          }
        }
        else {
          node->type = Node<Data>::Remote;
          node->n_children = BRANCH_FACTOR;
        }

        return false;
      }
      else {
        CkAbort("Error: can a light node be an internal node (above a TP root)?");
      }
    }
    */

    // create children
    node->n_children = BRANCH_FACTOR;
    Key child_key = node->key << LOG_BRANCH_FACTOR;
    int start = 0;
    int finish = start + node->n_particles;
    int non_local_children = 0;

    for (int i = 0; i < node->n_children; i++) {
      Key sibling_splitter = Utility::removeLeadingZeros(child_key + 1);

      // find number of particles in child
      int first_ge_idx;
      if (i < node->n_children - 1) {
        first_ge_idx = Utility::binarySearchGE(sibling_splitter, node->particles, start, finish);
      }
      else {
        first_ge_idx = finish;
      }
      int n_particles = first_ge_idx - start;

      /*
      // find owner treepieces of child
      int child_owner_start = owner_start;
      int child_owner_end;
      if (single_owner) {
        child_owner_end = child_owner_start;
      }
      else {
        if (i < node->n_children - 1) {
          int first_ge_tp = Utility::binarySearchGE(sibling_splitter, &splitters[0], owner_start, owner_end);
          child_owner_end = first_ge_tp - 1;
          owner_start = first_ge_tp;
        }
        else {
          child_owner_end = owner_end;
        }
      }
      */

      // create child and store in vector
      Node<Data>* child = new Node<Data>(child_key, node->depth + 1, n_particles, node->particles + start, 0, n_treepieces - 1, node);
      //Node<Data>* child = new Node<Data>(child_key, node->depth + 1, node->particles + start, n_particles, child_owner_start, child_owner_end, node);
      node->children.push_back(child);

      // recursive tree build
      bool local = recursiveBuild(child, saw_tp_key);

      if (!local) {
        non_local_children++;
      }

      start = first_ge_idx;
      child_key++;
    }
    if (non_local_children == 0) {
      node->type = Node<Data>::Internal;
    }
    else {
      node->type = Node<Data>::Boundary;
    }

    return (non_local_children == 0);
  }
  return false;
}

template <typename Data>
template <typename Visitor>
void TreePiece<Data>::upOnly(CProxy_TreeElement<Data> te_proxy) {
  global_data = te_proxy;

  std::queue<Node<Data>*> down, up;
  down.push(root);
  Visitor v (this->thisProxy, this->thisIndex);
  int n_curr_particles = 0;
  while (down.size()) {
    Node<Data>* node = down.front();
    down.pop();
    node->wait_count = 0;
    for (int i = 0; i < node->children.size(); i++) {
      switch (node->children[i]->type) {
        case Node<Data>::Internal:
        case Node<Data>::Boundary:
        case Node<Data>::Leaf: {
          down.push(node->children[i]);
          node->wait_count++;
        }
        default: break;
      }
    }
    if (node->type == Node<Data>::Leaf) {
      v.leaf(node);
      node->wait_count = 0;
      n_curr_particles += node->n_particles;
      up.push(node);
    }
  }
  if (!up.size()) { // no Data bc no Leaves, so send blank
    global_data[tp_key].template receiveData<Visitor>(this->thisProxy, Data(), this->thisIndex);
  }
  while (up.size()) {
    Node<Data>* node = up.front();
    up.pop();
    if (v.node(node)) { 
      if (--node->parent->wait_count == 0) up.push(node->parent);
    }
  }
}

template <typename Data>
template <typename Visitor>
void TreePiece<Data>::startDown(CProxy_CacheManager<Data> cache_manageri) {
  cache_manager = cache_manageri;
  root_from_tp_key = findNode(tp_key);
  if (cache_manager.ckLocalBranch()->connect(root_from_tp_key)) {
    cache_manager.ckLocalBranch()->template resumeTraversals<Visitor> (tp_key);
  }
  root_from_tp_key->parent->children[tp_key % 8] = NULL;
  root->triggerFree();

  curr_nodes = std::vector<std::set<Key> > (leaves.size());
  num_done = 0;
  trav_tops = std::vector<Node<Data>*> (leaves.size(), root);
  for (int i = 0; i < curr_nodes.size(); i++) {
    curr_nodes[i].insert(root->key);
  }
  goDown<Visitor>(root->key);
}

template <typename Data>
template <typename Visitor>
void TreePiece<Data>::startUpAndDown(CProxy_CacheManager<Data> cache_manageri) {
  cache_manager = cache_manageri;
  root_from_tp_key = findNode(tp_key);   
  if (cache_manager.ckLocalBranch()->connect(root_from_tp_key)) {
    cache_manager.ckLocalBranch()->template resumeTraversals<Visitor> (tp_key);
  }
  root_from_tp_key->parent->children[tp_key % 8] = NULL;
  root->triggerFree();

  if (!leaves.size()) {
    CkPrintf("tp %d finished!\n", this->thisIndex);
    CkCallback cb(CkReductionTarget(Main, doneDown), mainProxy);
    this->contribute(cb);
  }

  curr_nodes = std::vector<std::set<Key> > (leaves.size());
  num_done = 0;
  trav_tops = std::vector< Node<Data>* > (leaves.size(), NULL);
  for (int i = 0; i < curr_nodes.size(); i++) {
    curr_nodes[i].insert(leaves[i]->key);
    trav_tops[i] = leaves[i];
  }
  for (int i = 0; i < leaves.size(); i++) {
    goDown<Visitor> (leaves[i]->key);
  }
}

template<typename Data>
Node<Data>* TreePiece<Data>::findNode(Key key, bool start_from_tp_key) {
  std::vector<int> remainders;
  Key temp = key;
  while (temp >= 8 * (start_from_tp_key ? tp_key : 1)) {
    remainders.push_back(temp % 8);
    temp /= 8;
  }
  Node<Data>* node = start_from_tp_key ? root_from_tp_key : root;
  for (int i = remainders.size()-1; i >= 0; i--) {
    node = node->children[remainders[i]];
  }
  return node;
}

template <typename Data>
template <typename Visitor>
void TreePiece<Data>::requestNodes(Key key, CProxy_CacheManager<Data> cache_manager, int cm_index) {
  if (cm_index == CkMyPe()) return;

  Node<Data>* node = findNode(key, true);
  std::vector<Node<Data>> nodes;
  std::vector<Particle> sending_particles;
  node->tp_index = this->thisIndex;
  nodes.push_back(*node);
  for (int i = 0; i < node->n_children; i++) {
    Node<Data>* child = node->children[i];
    child->tp_index = this->thisIndex;
    nodes.push_back(*child);
    for (int j = 0; j < child->n_children; j++) {
      child->children[j]->tp_index = this->thisIndex;
      nodes.push_back(*(child->children[j]));
    }
  }
  for (int i = 0; i < nodes.size(); i++) {
    if (nodes[i].type == Node<Data>::Leaf) {
      for (int j = 0; j < nodes[i].n_particles; j++) {
        sending_particles.push_back(nodes[i].particles[j]);
      }
    }
  }
  MultiMsg<Data>* multimsg = new (sending_particles.size(), nodes.size()) MultiMsg<Data> (sending_particles.data(), sending_particles.size(), nodes.data(), nodes.size());
  cache_manager[cm_index].template addCache<Visitor>(multimsg);
}

template <typename Data>
template <typename Visitor>
void TreePiece<Data>::goDown(Key new_key) {
  Visitor v;
  if (new_key == 1) root = cache_manager.ckLocalBranch()->root;
  std::set<Key> to_go_down; 
  for (int i = 0; i < leaves.size(); i++) {
    if (curr_nodes[i].count(new_key) == 0) continue;
    curr_nodes[i].erase(new_key);
    std::stack<Node<Data>*> nodes;
    nodes.push(findNode(new_key));
    while (nodes.size()) {
      Node<Data>* node = nodes.top();
      nodes.pop();
      //CmiLock(cache_manager.ckLocalBranch()->clock);
      if (cache_manager.ckLocalBranch()->curr_waiting.count(node->key)) {
        //CkPrintf("non zero waiting count\n");
        curr_nodes[i].insert(node->key);
        if (cache_manager.ckLocalBranch()->curr_waiting[node->key].back() != this->thisIndex) // wont work on nodegroup
          cache_manager.ckLocalBranch()->curr_waiting[node->key].push_back(this->thisIndex);
        //CmiUnlock(cache_manager.ckLocalBranch()->clock);
        continue;
      }
      //CmiUnlock(cache_manager.ckLocalBranch()->clock);
      //CkPrintf("key = %d, type = %d\n", node->key, node->type);
      bool waiting = false;
      switch (node->type) {
        case Node<Data>::CachedBoundary:
        case Node<Data>::CachedRemote:
        case Node<Data>::Internal: {
          if (v.node(node, leaves[i])) {
            for (int i = 0; i < node->children.size(); i++) {
              nodes.push(node->children[i]);
            }
          }
          break;
        }
        case Node<Data>::CachedRemoteLeaf:
        case Node<Data>::Leaf: {
          v.leaf(node, leaves[i]);
          break;
        } 
        case Node<Data>::Boundary:
        case Node<Data>::RemoteAboveTPKey: {
          global_data[node->key].template requestData<Visitor>(cache_manager, CkMyPe());
          waiting = true;
          break;
        }
        case Node<Data>::Remote:
        case Node<Data>::RemoteLeaf: {
          this->thisProxy[node->tp_index].template requestNodes<Visitor>(node->key, cache_manager, CkMyPe());
          waiting = true;          
        }
        default: break; 
      }
      if (waiting) {
        curr_nodes[i].insert(node->key);
        std::vector<int> myvec (1, this->thisIndex);
        //CmiLock(cache_manager.ckLocalBranch()->clock);
        cache_manager.ckLocalBranch()->curr_waiting.insert(make_pair(node->key, myvec));
        //CmiUnlock(cache_manager.ckLocalBranch()->clock);
      }
    } 
    if (!curr_nodes[i].size()) {
      if (trav_tops[i]->parent == NULL) {
        num_done++;
      }
      else {        
        for (int j = 0; j < trav_tops[i]->parent->children.size(); j++) {
          Node<Data>* child = trav_tops[i]->parent->children[j];
          if (child != trav_tops[i]) {
             if (trav_tops[i]->parent->type == Node<Data>::Boundary) child->type = Node<Data>::RemoteAboveTPKey;
             curr_nodes[i].insert(child->key);
             to_go_down.insert(child->key);
          }
        }
        trav_tops[i] = trav_tops[i]->parent;
      }
    }
  }
  if (num_done == leaves.size()) {
    CkPrintf("tp %d finished!\n", this->thisIndex);
    CkCallback cb(CkReductionTarget(Main, doneDown), mainProxy);
    this->contribute(cb);
  }
  for (std::set<Key>::iterator it = to_go_down.begin(); it != to_go_down.end(); it++) {
    goDown<Visitor> (*it);
  } 
}

template <typename Data>
void TreePiece<Data>::perturb (Real timestep) {
  for (int i = 0; i < leaves.size(); i++) {
    for (int i = 0; i < leaves[i]->n_particles; i++) {
      leaves[i]->particles[i].perturb(timestep, leaves[i]->data.sum_forces[i]);
    }
  }
}

template <typename Data>
void TreePiece<Data>::flush(CProxy_Reader readers) {
  // debug
  flushed_particles.resize(0);
  flushed_particles.insert(flushed_particles.end(), particles.begin(), particles.end());

  ParticleMsg *msg = new (particles.size()) ParticleMsg(particles.data(), particles.size());
  readers[CkMyPe()].receive(msg);
  particles.resize(0);
  particle_index = 0;
}

template <typename Data>
void TreePiece<Data>::rebuild(const CkCallback& cb) {
  // clean up old tree information
  root->triggerFree();
  root = nullptr;

  // build a new tree from scratch
  build(cb);
}

template <typename Data>
void TreePiece<Data>::print(Node<Data>* root) {
  ostringstream oss;
  oss << "tree." << this->thisIndex << ".dot";
  ofstream out(oss.str().c_str());
  CkAssert(out.is_open());
  out << "digraph tree" << this->thisIndex << "{" << endl;
  root->dot(out);
  out << "}" << endl;
  out.close();
}

#endif // SIMPLE_TREEPIECE_H_
