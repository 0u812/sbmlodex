import tesbml as libsbml

import os

def getRoadrunnerSBML(obj):
    """Checks if obj is a Roadrunner instance. Returns the SBML if true, else raises TypeError."""
    try:
        from roadrunner import RoadRunner
    except ImportError:
        raise TypeError('Not a Roadrunner instance')

    if isinstance(obj, RoadRunner):
        return obj.getSBML()
    else:
        raise TypeError('Not a Roadrunner instance')

def readSBMLModel(sbml_model):
    try:
        sbml = getRoadrunnerSBML(sbml_model)
    except TypeError:
        if os.path.exists(sbml_model) and os.path.isfile(sbml_model):
            # it's a file path
            if os.path.splitext(sbml_model)[1] == '.sb':
                # it's an Antimony file
                with open(sbml_model) as f:
                    return libsbml.readSBML(antimonyConverter.antimonyToSBML(f.read()))
            elif os.path.splitext(sbml_model)[1] == '.txt':
                raise RuntimeError('File ending in ".txt" is ambiguous - pass an SBML file (.xml) or an Antimony file (.sb).')
            else:
                return libsbml.readSBMLFromFile(sbml_model)
        else:
            # check if it's Antimony source
            try:
                return libsbml.readSBML(antimonyConverter.antimonyToSBML(sbml_model))
            except:
                # it better be SBML
                # this will throw if it's not SBML
                return libsbml.readSBML(model)

def extractODEsToString(sbml_model):
    """Extracts a list of ODEs from sbml_model. sbml_model can be either a Roadrunner
    instance, a path to an SBML file, an SBML string, or an Antimony string."""
    extractor = ODEExtractor(readSBMLModel(sbml_model))
    return extractor.toString()

class Accumulator:
    def __init__(self, species_id):
        self.reaction_map = {}
        self.reactions = []
        self.species_id = species_id

    def addReaction(self, reaction, stoich):
        rid = reaction.getId()
        if rid in self.reaction_map:
            self.reaction_map[rid]['stoich'] += stoich
        else:
            self.reaction_map[rid] = {
                'reaction': reaction,
                'formula': self.getFormula(reaction),
                'stoich': stoich,
            }
            self.reactions.append(rid)

    def getFormula(self, reaction):
        return reaction.getKineticLaw().getFormula()

    def toString(self):
        lhs = 'd{}/dt'.format(self.species_id)
        terms = []
        for rid in self.reactions:
            if self.reaction_map[rid]['stoich'] != 1:
                stoich = ''
            else:
                stoich = str(abs(self.reaction_map[rid]['stoich'])) + '*'

            if len(terms) > 0:
                if self.reaction_map[rid]['stoich'] < 0:
                    op = ' - '
                else:
                    op = ' + '
            else:
                if self.reaction_map[rid]['stoich'] < 0:
                    op = '-'
                else:
                    op = ''
            terms.append(op + stoich + self.reaction_map[rid]['formula'])
        rhs = ''.join(terms)
        return lhs + ' = ' + rhs

class ODEExtractor:
    def __init__(self, sbml_doc):
        self.sbml_doc = sbml_doc

        if self.sbml_doc.getNumErrors() > 0:
            print('{} errors while reading SBML'.format(self.sbml_doc.getNumErrors()))
            for e in (self.sbml_doc.getError(i) for i in range(self.sbml_doc.getNumErrors())):
                print(e.getMessage())
            raise RuntimeError('Errors reading SBML')

        self.model = sbml_doc.getModel()

        self.species_map = {}
        self.species_symbol_map = {}
        self.use_species_names = False

        from collections import defaultdict
        self.accumulators = {}
        self.accumulator_list = []

        def reactionParticipant(participant, stoich):
            stoich_sign = 1
            if stoich < 0:
                stoich_sign = -1
            if participant.isSetStoichiometry():
                stoich = participant.getStoichiometry()
            elif isSetStoichiometryMath():
                raise RuntimeError('Stoichiometry math not supported')
            self.accumulators[participant.getSpecies()].addReaction(r, stoich_sign*stoich)

        newReactant = lambda p: reactionParticipant(p, 1)

        newProduct  = lambda p: reactionParticipant(p, -1)

        for s in (self.model.getSpecies(i) for i in range(self.model.getNumSpecies())):
            self.species_map[s.getId()] = s
            if s.isSetName() and self.use_species_names:
                self.species_symbol_map[s.getId()] = s.getName()
            else:
                self.species_symbol_map[s.getId()] = s.getId()
            a = Accumulator(s.getId())
            self.accumulators[s.getId()] = a
            self.accumulator_list.append(a)

        for r in (self.model.getReaction(i) for i in range(self.model.getNumReactions())):
            for reactant in (r.getReactant(i) for i in range(r.getNumReactants())):
                newReactant(reactant)
            for product in (r.getProduct(i) for i in range(r.getNumProducts())):
                newProduct(product)

    def toString(self):
        r = ''
        for a in self.accumulator_list:
            r += a.toString() + '\n'
        return r
