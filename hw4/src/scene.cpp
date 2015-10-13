/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2012 by Wenzel Jakob and Steve Marschner.

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/luminaire.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) 
	: m_integrator(NULL), m_sampler(NULL), m_camera(NULL), 
	  m_medium(NULL), m_envLuminaire(NULL), m_evaluator(NULL) {
	m_kdtree = new KDTree();
}

Scene::~Scene() {
	delete m_kdtree;
	if (m_sampler)
		delete m_sampler;
	if (m_camera)
		delete m_camera;
	if (m_integrator)
		delete m_integrator;
	if (m_medium)
		delete m_medium;
	if (m_envLuminaire)
		delete m_envLuminaire;
}

Color3f Scene::sampleDirect(LuminaireQueryRecord &lRec, const Point2f &_sample) const {
	if (m_luminaires.size() == 0)
		throw NoriException("Scene::sampleDirect(): No luminaires were defined!");

	Point2f sample(_sample);
	int index = std::min((int) (m_luminaires.size() * sample.x()), (int) m_luminaires.size()-1);
	sample.x() = m_luminaires.size() * sample.x() - index;

	lRec.luminaire = m_luminaires[index];
	Color3f value = lRec.luminaire->sample(lRec, sample);

	if (lRec.pdf != 0) {
		if (rayIntersect(Ray3f(lRec.ref, lRec.d, Epsilon, lRec.dist * (1-1e-4f))))
			return Color3f(0.0f);
		lRec.pdf /= m_luminaires.size();
		return value * m_luminaires.size();
	} else {
		return Color3f(0.0f);
	}
}

float Scene::pdfDirect(const LuminaireQueryRecord &lRec) const {
	return lRec.luminaire->pdf(lRec) / m_luminaires.size();
}

bool Scene::sampleDistance(const Ray3f &ray, Sampler *sampler, float &t, Color3f &weight) const {
	if (m_medium) {
		return m_medium->sampleDistance(ray, sampler, t, weight);
	} else {
		weight = Color3f(1.0f);
		return false;
	}
}

Color3f Scene::evalTransmittance(const Ray3f &ray, Sampler *sampler) const {
	if (m_medium) {
		return m_medium->evalTransmittance(ray, sampler);
	} else {
		return Color3f(1.0f);
	}
}

void Scene::activate() {
	m_kdtree->build();

	if (!m_integrator)
		throw NoriException("No integrator was specified!");
	if (!m_camera)
		throw NoriException("No camera was specified!");
	
	if (!m_sampler) {
		/* Create a default (independent) sampler */
		m_sampler = static_cast<Sampler*>(
			NoriObjectFactory::createInstance("independent", PropertyList()));
	}

	cout << endl;
	cout << "Configuration: " << qPrintable(toString()) << endl;
	cout << endl;
}

void Scene::addChild(NoriObject *obj) {
	switch (obj->getClassType()) {
		case EMesh: {
				Mesh *mesh = static_cast<Mesh *>(obj);
				m_kdtree->addMesh(mesh);
				m_meshes.push_back(mesh);
				if (mesh->isLuminaire())
					m_luminaires.push_back(mesh->getLuminaire());
			}
			break;
		
		case ELuminaire: {
				Luminaire *luminaire = static_cast<Luminaire *>(obj);
				if (!luminaire->isEnvironmentLuminaire())
					throw NoriException("The specified root-level luminaire is not an environment luminaire!");
				if (m_envLuminaire)
					throw NoriException("There can only be one environment luminaire per scene!");
				m_envLuminaire = luminaire;
				m_luminaires.push_back(luminaire);
			}
			break;
                        
                case EEvaluator: {
				Evaluator *eval = static_cast<Evaluator *>(obj);
				if (m_envLuminaire)
					throw NoriException("There can only be one scene evaluator!");
				m_evaluator = eval;
			}
			break;  

		case ESampler:
			if (m_sampler)
				throw NoriException("There can only be one sampler per scene!");
			m_sampler = static_cast<Sampler *>(obj);
			break;

		case ECamera:
			if (m_camera)
				throw NoriException("There can only be one camera per scene!");
			m_camera = static_cast<Camera *>(obj);
			break;
		
		case EMedium:
			if (m_medium)
				throw NoriException("There can only be one medium per scene!");
			m_medium = static_cast<Medium *>(obj);
			break;

		case EIntegrator:
			if (m_integrator)
				throw NoriException("There can only be one integrator per scene!");
			m_integrator = static_cast<Integrator *>(obj);
			break;

		default:
			throw NoriException(QString("Scene::addChild(<%1>) is not supported!").arg(
				classTypeName(obj->getClassType())));
	}
}

QString Scene::toString() const {
	QString meshes;
	for (size_t i=0; i<m_meshes.size(); ++i) {
		meshes += QString("  ") + indent(m_meshes[i]->toString(), 2);
		if (i + 1 < m_meshes.size())
			meshes += ",";
		meshes += "\n";
	}
	return QString(
		"Scene[\n"
		"  integrator = %1,\n"
		"  sampler = %2\n"
		"  camera = %3,\n"
		"  medium = %4,\n"
		"  envLuminaire = %5,\n"
		"  meshes = {\n"
		"  %6},\n"
                "  evaluator = %7\n"
		"]")
	.arg(indent(m_integrator->toString()))
	.arg(indent(m_sampler->toString()))
	.arg(indent(m_camera->toString()))
	.arg(m_medium ? indent(m_medium->toString()) : QString("null"))
	.arg(indent(m_envLuminaire ? m_envLuminaire->toString() : QString("null")))
	.arg(indent(meshes, 2))
        .arg(indent(m_evaluator ? m_evaluator->toString() : QString("null")));
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
