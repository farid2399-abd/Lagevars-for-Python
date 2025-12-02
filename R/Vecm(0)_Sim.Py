# simulation.py
import numpy as np
from typing import Tuple, Optional

def simular_vecm_configurable(
    T: int = 200,
    N: int = 10,
    r: int = 1,
    p: int = 0,
    seed: Optional[int] = None,
    alpha_true: Optional[np.ndarray] = None,
    beta_true: Optional[np.ndarray] = None,
    sigma_epsilon: float = 0.2
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simula un proceso VECM(p) con parámetros configurables
    
    Parámetros:
    -----------
    T : int
        Número total de observaciones temporales
    N : int
        Número de series temporales
    r : int
        Rango de cointegración (número de relaciones de cointegración)
    p : int
        Orden del VECM (0 para VECM(0))
    seed : int, opcional
        Semilla para reproducibilidad
    alpha_true : np.ndarray, opcional
        Matriz alpha verdadera de forma (N, r)
    beta_true : np.ndarray, opcional
        Matriz beta verdadera de forma (N, r)
    sigma_epsilon : float
        Desviación estándar de los errores
    
    Retorna:
    --------
    X : np.ndarray
        Matriz de datos simulados (T, N)
    beta_true : np.ndarray
        Vector de cointegración verdadero
    alpha_true : np.ndarray
        Matriz de ajuste verdadera
    Pi : np.ndarray
        Matriz Pi = alpha * beta'
    
    Ejemplo:
    --------
    >>> X, beta, alpha, Pi = simular_vecm_configurable(T=200, N=10, r=1, p=0)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Si no se proporcionan parámetros verdaderos, generarlos
    if beta_true is None:
        # Solo las primeras r series están cointegradas
        beta_true = np.zeros(N)
        if r > 0:
            beta_true[:r] = np.random.choice([1, -1], r)
            # Normalizar: primer elemento = 1
            if beta_true[0] != 0:
                beta_true = beta_true / beta_true[0]
    
    if alpha_true is None:
        alpha_true = np.zeros(N)
        if r > 0:
            # Valores típicos para alpha (entre -0.8 y 0.8)
            alpha_true[:r] = np.random.uniform(-0.8, 0.8, r)
    
    # Asegurar que son vectores
    beta_true = np.asarray(beta_true).flatten()
    alpha_true = np.asarray(alpha_true).flatten()
    
    # Matriz Pi = alpha * beta'
    Pi = np.outer(alpha_true, beta_true)
    
    # Simular serie según VECM(p)
    X = np.zeros((T, N))
    
    # Valores iniciales
    X[0, :] = np.random.normal(0, 1, N)
    
    # Simular errores
    epsilon = np.random.normal(0, sigma_epsilon, (T, N))
    
    # Coeficientes para términos de diferencia retardados (si p > 0)
    if p > 0:
        # Generar coeficientes pequeños para los términos de diferencia
        gamma = np.random.uniform(-0.2, 0.2, (p, N, N)) / (np.arange(1, p+1)[:, np.newaxis, np.newaxis])
    else:
        gamma = None
    
    # Simular proceso VECM
    for t in range(1, T):
        # Término de corrección de error
        correction = Pi @ X[t-1, :]
        
        # Agregar términos de diferencia retardados si p > 0
        if p > 0:
            for lag in range(1, min(p, t) + 1):
                correction += gamma[lag-1] @ (X[t-lag, :] - X[t-lag-1, :] if t > lag else X[t-lag, :])
        
        # Actualizar
        X[t, :] = X[t-1, :] + correction + epsilon[t, :]
    
    return X, beta_true, alpha_true, Pi


def simular_vecm_ec3(T: int, N: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Versión específica para la ecuación 3 del artículo (VECM(0) con estructura conocida)
    
    Parámetros:
    -----------
    T : int
        Número de observaciones
    N : int
        Número de series
    
    Retorna:
    --------
    X : np.ndarray
        Datos simulados
    beta_true : np.ndarray
        Beta verdadero [1, -1, 0, ..., 0]
    alpha_true : np.ndarray
        Alpha verdadero [-0.5, 0.3, 0, ..., 0]
    Pi : np.ndarray
        Matriz Pi = alpha * beta'
    """
    # Parámetros fijos según el artículo
    beta_true = np.array([1, -1] + [0]*(N-2))
    alpha_true = np.array([-0.5, 0.3] + [0]*(N-2))
    
    return simular_vecm_configurable(
        T=T, N=N, r=1, p=0,
        beta_true=beta_true,
        alpha_true=alpha_true,
        sigma_epsilon=0.2
    )
