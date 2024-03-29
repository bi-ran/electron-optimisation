#ifndef COMMON_H
#define COMMON_H

#define DECLARE(type, var) type var;
#define DECLPTR(type, var) type* var;

#define ZERO(type, var) var = 0;
#define INVALID(type, var) var = -1;
#define NEWVEC(type, var) var = new type();

#define CLEAR(type, var) var->clear();

#define RREF(type, var, tree)                                               \
    tree->SetBranchStatus(#var, 1);                                         \
    tree->SetBranchAddress(#var, &var);
#define RVAR(type, var, tree)                                               \
    tree->SetBranchStatus(#var, 1);                                         \
    tree->SetBranchAddress(#var, var);

#define BRNREF(type, var, tree) tree->Branch(#var, &var);
#define BRNVAR(type, var, tree) tree->Branch(#var, var);

#define VARCOPY(type, var, tree)                                            \
    var = tree->var;

#define VECCOPY(type, var, tree)                                            \
    if (tree->var != nullptr)                                               \
        std::copy(tree->var->begin(), tree->var->end(),                     \
                  std::back_inserter(*var));

#endif /* COMMON_H */
