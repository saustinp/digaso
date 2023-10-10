#ifndef __DMD_H
#define __DMD_H

void createdmdstruct(dmdstruct &dmd, meshstruct &mesh);

void setupsys(sysstruct &sys, meshstruct &mesh, appstruct &app);

void writedmdstruct(string filename, dmdstruct &dmd);

#endif
