#include <PeleLMDerive.H>

using namespace amrex;

PeleLMDeriveRec::PeleLMDeriveRec (const std::string& a_name,
                                  IndexType          result_type,
                                  int                nvar_derive,
                                  PeleLMDeriveFunc   der_func,
                                  DeriveBoxMap       box_map,
                                  Interpolater*      a_interp)
    :
    derive_name(a_name),
    variable_names(),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(box_map)
{}

PeleLMDeriveRec::PeleLMDeriveRec (const std::string&   a_name,
                                  IndexType            result_type,
                                  int                  nvar_derive,
		                            Vector<std::string>& var_names,
                                  PeleLMDeriveFunc     der_func,
                                  DeriveBoxMap         box_map,
                                  Interpolater*        a_interp)
    :
    derive_name(a_name),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(box_map)
{}

PeleLMDeriveRec::~PeleLMDeriveRec () 
{
   func     = nullptr;
   mapper   = 0;
   bx_map   = 0;
}

const std::string&
PeleLMDeriveRec::name () const noexcept
{
    return derive_name;
}

IndexType
PeleLMDeriveRec::deriveType () const noexcept
{
    return der_type;
}

PeleLMDeriveFunc
PeleLMDeriveRec::derFunc () const noexcept
{
    return func;
}

PeleLMDeriveRec::DeriveBoxMap
PeleLMDeriveRec::boxMap () const noexcept
{
    return bx_map;
}

Interpolater*
PeleLMDeriveRec::interp () const noexcept
{
    return mapper;
}

int
PeleLMDeriveRec::numDerive () const noexcept
{
    return n_derive;
}

const
std::string&
PeleLMDeriveRec::variableName(int comp) const noexcept
{
  if (comp < variable_names.size()) 
     return variable_names[comp];

  return derive_name;
}

PeleLMDeriveList::PeleLMDeriveList () {}

void
PeleLMDeriveList::add (const std::string&            name,
                       IndexType                     result_type,
                       int                           nvar_der,
                       PeleLMDeriveFunc              der_func,
                       PeleLMDeriveRec::DeriveBoxMap bx_map,
                       Interpolater*                 interp)
{
    lst.push_back(PeleLMDeriveRec(name,result_type,nvar_der,der_func,bx_map,interp));
}

void
PeleLMDeriveList::add (const std::string&            name,
                       IndexType                     result_type,
                       int                           nvar_der,
                       Vector<std::string>&          vars,
                       PeleLMDeriveFunc              der_func,
                       PeleLMDeriveRec::DeriveBoxMap bx_map,
                       Interpolater*                 interp)
{
    lst.push_back(PeleLMDeriveRec(name,result_type,nvar_der,vars,der_func,bx_map,interp));
}

std::list<PeleLMDeriveRec>&
PeleLMDeriveList::dlist ()
{
    return lst;
}

bool
PeleLMDeriveList::canDerive (const std::string& name) const 
{
    for (std::list<PeleLMDeriveRec>::const_iterator li = lst.begin(), End = lst.end();
         li != End;
         ++li)
    {
        if (li->derive_name == name)
            return true;
    }
    return false;
}

const PeleLMDeriveRec*
PeleLMDeriveList::get (const std::string& name) const
{
    for (std::list<PeleLMDeriveRec>::const_iterator li = lst.begin(), End = lst.end();
         li != End;
         ++li)
    {
        if (li->derive_name == name)
            return &(*li);
    }
    return 0;
}
