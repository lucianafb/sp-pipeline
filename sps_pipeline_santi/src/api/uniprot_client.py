"""
UniProt REST API Client con paginación completa.

Maneja cursor-based pagination, retry con backoff, y barra de progreso.
"""

import re
import time
import logging
import requests
from typing import Generator, Optional

from tqdm import tqdm

logger = logging.getLogger(__name__)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# Timeout por request (connect, read) en segundos
DEFAULT_TIMEOUT = (10, 60)
MAX_RETRIES = 3
BACKOFF_FACTOR = 2  # 2s, 4s, 8s


def _parse_next_link(link_header: str) -> Optional[str]:
    """Extrae la URL 'next' del header Link de UniProt."""
    if not link_header:
        return None
    # Formato: <URL>; rel="next"
    match = re.search(r'<([^>]+)>;\s*rel="next"', link_header)
    return match.group(1) if match else None


def _request_with_retry(url: str, params: Optional[dict] = None) -> requests.Response:
    """Hace un GET con retry y backoff exponencial en errores 429/5xx."""
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, params=params, timeout=DEFAULT_TIMEOUT)

            if response.status_code == 200:
                return response

            if response.status_code == 429:
                # Rate limited — respetar Retry-After si viene
                retry_after = int(response.headers.get("Retry-After", BACKOFF_FACTOR ** (attempt + 1)))
                logger.warning(f"Rate limited (429). Esperando {retry_after}s...")
                time.sleep(retry_after)
                continue

            if response.status_code >= 500:
                wait = BACKOFF_FACTOR ** (attempt + 1)
                logger.warning(f"Error del servidor ({response.status_code}). Retry en {wait}s...")
                time.sleep(wait)
                continue

            # 4xx que no es 429: error del cliente, no reintentamos
            response.raise_for_status()

        except requests.exceptions.ConnectionError as e:
            wait = BACKOFF_FACTOR ** (attempt + 1)
            logger.warning(f"Error de conexión: {e}. Retry en {wait}s...")
            time.sleep(wait)

        except requests.exceptions.Timeout:
            wait = BACKOFF_FACTOR ** (attempt + 1)
            logger.warning(f"Timeout. Retry en {wait}s...")
            time.sleep(wait)

    raise RuntimeError(f"Falló después de {MAX_RETRIES} intentos: {url}")


def fetch_entries(
    query: str,
    max_results: Optional[int] = None,
    page_size: int = 500,
) -> Generator[dict, None, None]:
    """
    Genera entradas de UniProt paginando automáticamente.

    Args:
        query: Query string para la API de UniProt.
        max_results: Límite opcional de resultados (None = todos).
        page_size: Resultados por página (max 500).

    Yields:
        dict: Cada entrada JSON de UniProt tal cual viene de la API.
    """
    page_size = min(page_size, 500)  # UniProt cap

    params = {
        "query": query,
        "format": "json",
        "size": page_size,
    }

    logger.info(f"Iniciando query: {query}")

    # Primera request
    response = _request_with_retry(UNIPROT_SEARCH_URL, params=params)
    data = response.json()

    total_results = int(response.headers.get("x-total-results", 0))
    effective_total = total_results
    if max_results is not None:
        effective_total = min(total_results, max_results)

    logger.info(f"Total en UniProt: {total_results} | A recuperar: {effective_total}")

    yielded = 0

    with tqdm(total=effective_total, desc="Descargando", unit=" entries") as pbar:
        while True:
            hits = data.get("results", [])

            for entry in hits:
                if max_results is not None and yielded >= max_results:
                    return
                yield entry
                yielded += 1
                pbar.update(1)

            # Buscar siguiente página
            link_header = response.headers.get("Link", "")
            next_url = _parse_next_link(link_header)

            if not next_url:
                break  # No hay más páginas

            if max_results is not None and yielded >= max_results:
                break

            # Siguiente página (URL completa, sin params adicionales)
            response = _request_with_retry(next_url)
            data = response.json()

    logger.info(f"Descarga completa: {yielded} entradas recuperadas.")


def fetch_all_entries(
    query: str,
    max_results: Optional[int] = None,
) -> list[dict]:
    """
    Función de conveniencia que retorna todas las entradas como lista.

    Para datasets grandes, preferir fetch_entries() como generador.
    """
    return list(fetch_entries(query, max_results=max_results))
